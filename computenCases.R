
#' Sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement
#' error and treatment efficacy
#'
#' Performs sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement error and treatment efficacy [Gilbert, Janes, and Huang (2015).
#' ``Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials.'']
#' 
#' @param tau Month at which immune responses are measured (beginning of follow-up for the endpoint)
#' @param taumax Month at which follow-up for the endpoint ends
#' @param VEoverall Overal cumulative VE between month tau and month taumax
#' @param NvEnrolled Number enrolled in the subunit vaccine group
#' @param N Number of subjects in the vaccine group at risk at month tau
#' @param risk0 Cumulative failure rate in the hypothetical placebo arm between Month tau and taumax
#' @param dropOutRate Cumulative dropout rate between month 0 and taumax
#' @param fracImmunog Fraction of enrolled cases selected for the immunogenicity subset
#' 
#' @details
#' This function performs simulations to calculate the power for testing whether a trichotomous or continuous biomarker
#' is a correlate of risk (CoR).
#' 
#' N can be computed using the function `computeN`
#' 
#'   Calculate the number of observed primary endpoint cases in a vaccine group 
#' registered between Month tau and taumax, which is called n_cases in the paper.  
#' Calculate this under the following assumptions:
#'
#' 1. The failure time T and the censoring time C are independent.
#' 2. T and C have exponential distributions in the hypothetical placebo arm
#'    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
#' 3. RRoverall = P(tau < T <= t|Z=1)/P(tau < T <= t|Z=0) for all t between
#'                tau and taumax (this will only approximately hold).
#'
#' Number of observed rotavirus endpoints between Month tau and taumax in vaccinees:  
#' Logic: = Number vaccinees at-risk at Mo tau * P(T<=taumax,T<=C|X>tau,Z=1)
#'        = \eqn{N*\int_0^infty P(T<=min(c,taumax),C=c|X>tau,Z=1)dc}
#'        = \eqn{N*\int_tau^infty P(tau<T<=min(c,taumax),C=c|Z=1)dc / P(X>tau|Z=1)}
#'        = \eqn{N*{\int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc}/P(X>tau|Z=1)}.
#'
#' Now, the above integral in { } can be written as
#'   \eqn{\int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc}
#' = \eqn{\int_tau^taumax P(tau<T<=c|Z=1)P(C=c)dc + P(tau<T<=taumax|Z=1)\int_taumax^infty P(C=c)dc}
#' = \eqn{RRoverall*\int_tau^taumax{[exp(-thetat*tau)-exp(-thetat*c)]
#'            *thetac*exp(-thetac*c)dc} (term1)
#' + RRoverall*[pexp(taumax,thetat)-pexp(tau,thetat)]*[1 - pexp(taumax,thetac)].} (term2)
#'
#' Term (1) is calculated to be:
#'   \eqn{RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
#' - thetac*\int_tau^taumax exp(-(thetat+thetac)c)dc]}
#' = \eqn{RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
#' - (thetac/(thetac+thetat))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))}
#'
#' Lastly, the denominator (call it term3) equals:
#' P(X>tau|Z=1) = \eqn{[1-RRoverall0toTau*P(T<=tau|Z=0)]P(C>tau)}
#'              = \eqn{[1-RRoverall0toTau*pexp(tau,thetat)][1-pexp(tau,thetac)].}
#'
#' @return List with the following elements:
#'   nAtRiskTauCases (nCases in manuscript): Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax.
#    nAtRiskTauControls (nControls in manusript): Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
#    nAtRiskTauCasesPhase2 (nCasesPhase2): Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax and with the biomarker measured
#'
#' @examples
#'
#' @import survival
#' @import osDesign
#' @export
computenCases <- function(tau, taumax, VEoverall, NvEnrolled, N, risk0, dropOutRate, fracImmunog) {

  # With T the failure time, C the censoring time, X = min(T,C), and 
  # Z vaccination status (vaccine or hypothetical placebo), we have
  # RRoverall = P(tau < T <= taumax|Z=1)/P(tau < T <= taumax|Z=0) 
  RRoverall <- 1 - VEoverall

  # Cumulative failure rate in the hypothetical placebo arm
  # between Month 0 and taumax:
  risk0Mo0toTaumax <- risk0*(1 + tau/taumax)
  
  # Specify VE in the early period between enrollment and time tau --
  # specify it half way between 0 and VEoverall:
  VEoverall0toTau <- VEoverall/2
  RRoverall0toTau <- 1 - VEoverall0toTau  
  
  # Calculate the number of observed primary endpoint cases in a vaccine group 
  # registered between Month tau and taumax, which is called n_cases in the paper.  
  # Calculate this under the following assumptions:
  #
  # 1. The failure time T and the censoring time C are independent.
  # 2. T and C have exponential distributions in the hypothetical placebo arm
  #    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
  # 3. RRoverall = P(tau < T <= t|Z=1)/P(tau < T <= t|Z=0) for all t between
  #                tau and taumax (this will only approximately hold).
  #
  
  # Calculate thetat and thetac on a monthly time-scale:
  thetat <- -log(1-risk0Mo0toTaumax)/taumax
  thetac <- -log(1-dropOutRate)/taumax
  
  # Number of observed rotavirus endpoints between Month tau and taumax in vaccinees:  
  # Logic: = Number vaccinees at-risk at Mo tau * P(T<=taumax,T<=C|X>tau,Z=1)
  #        = N*\int_0^infty P(T<=min(c,taumax),C=c|X>tau,Z=1)dc
  #        = N*\int_tau^infty P(tau<T<=min(c,taumax),C=c|Z=1)dc / P(X>tau|Z=1)
  #        = N*{\int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc}/P(X>tau|Z=1).
  
  # Now, the above integral in { } can be written as
  #   \int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc
  # = \int_tau^taumax P(tau<T<=c|Z=1)P(C=c)dc + P(tau<T<=taumax|Z=1)\int_taumax^infty P(C=c)dc
  # = RRoverall*\int_tau^taumax{[exp(-thetat*tau)-exp(-thetat*c)]
  #            *thetac*exp(-thetac*c)dc} (term1)
  # + RRoverall*[pexp(taumax,thetat)-pexp(tau,thetat)]*[1 - pexp(taumax,thetac)]. (term2)
  
  # Term (1) is calculated to be:
  #   RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
  # - thetac*\int_tau^taumax exp(-(thetat+thetac)c)dc]
  # = RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
  # - (thetac/(thetac+thetat))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))
  
  # Lastly, the denominator (call it term3) equals:
  # P(X>tau|Z=1) = [1-RRoverall0toTau*P(T<=tau|Z=0)]P(C>tau)
  #              = [1-RRoverall0toTau*pexp(tau,thetat)][1-pexp(tau,thetac)].
  
  # Putting it together
  term1 <- RRoverall*(exp(-thetat*tau)*(pexp(taumax,thetac)-pexp(tau,thetac))
                      - (thetac/(thetac+thetat))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))
  
  term2 <- RRoverall*(pexp(taumax,thetat)-pexp(tau,thetat))*(1-pexp(taumax,thetac))
  term3 <- (1-RRoverall0toTau*pexp(tau,thetat))*(1-pexp(tau,thetac))
  
  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax.
  # nCases is notation in the paper. nAtRiskTauCases is notation used in R program computepower
  nAtRiskTauCases <- round(N*(term1+term2)/term3)
  
  # Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
  # Calculated the same as for term3 except using taumax instead of tau
  # nControls is notation in the paper. nAtRiskTauControls is notation used in R program computepower
  nAtRiskTauControls <- round(NvEnrolled*(1-RRoverall*pexp(taumax,thetat))*(1-pexp(taumax,thetac)))
  
  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax 
  # and with the biomarker measured (i.e., in Phase 2).
  # nAtRiskTauCasesPhase2 is notation used in R program computepower
  nAtRiskTauCasesPhase2 <- round(fracImmunog*nAtRiskTauCases)
  
  return(list(nAtRiskTauCases = nAtRiskTauCases, nAtRiskTauControls = nAtRiskTauControls, nAtRiskTauCasesPhase2 = nAtRiskTauCasesPhase2))
}

