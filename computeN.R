#' Sample size calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement
#' error and treatment efficacy
#'
#' Performs sample size/power calculations for assessing biomarkers as correlates of risk (CoRs) accounting for measurement error and treatment efficacy [Gilbert, Janes, and Huang (2015).
#'
#' @param tau Month at which immune responses are measured (beginning of follow-up for the endpoint)
#' @param taumax Month at which follow-up for the endpoint ends
#' @param VEoverall Overal cumulative VE between month tau and month taumax
#' @param NvEnrolled Number enrolled in the subunit vaccine group
#' @param risk0 Cumulative failure rate in the hypothetical placebo arm between Month tau and taumax
#' @param dropOutRate Cumulative dropout rate between month 0 and taumax
#' @param fracImmunog Fraction of enrolled cases selected for the immunogenicity subset
#' 
#' @details
#' This function performs simulations to calculate the power for testing whether a trichotomous or continuous biomarker
#' is a correlate of risk (CoR).
#'
#' @return List with the following elements: 
#'   N - the number of subjects in the vaccine group at risk at month tau
#'   nCases (nCases in manuscript): Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax.
#    nControls (nControls in manusript): Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
#    nCasesWithS (nCasesPhase2): Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax and with the biomarker measured
#'     
#'
#' @examples
#' tau <- 3.5
#' taumax <- 24
#' VEoverall <- 0.75
#' NvEnrolled <- 4100
#' risk0 <- 0.034
#' dropOutRate <- 0.1
#' fracImmunog <- 1.0
#' computeN(tau, taumax, VEoverall, NvEnrolled, risk0, dropOutRate, fracImmunog)
#' 
#'
#' @import survival
#' @import osDesign
#' @export
computeN <- function(tau, taumax, VEoverall, NvEnrolled, risk0, dropOutRate, fracImmunog) {
  
  # With T the failure time, C the censoring time, X = min(T,C), and 
  # Z vaccination status (vaccine or hypothetical placebo), we have
  # RRoverall = P(tau < T <= taumax|Z=1)/P(tau < T <= taumax|Z=0) 
  RRoverall <- 1 - VEoverall
  
  # Cumulative failure rate in the hypothetical placebo arm between Month 0 and 24:
  risk0Mo0toTaumax <- risk0*(1 + tau/taumax)
  
  # Specify VE in the early period between enrollment and time tau 
  # Specify it half way between 0 and VEoverall:
  VEoverall0toTau <- VEoverall/2
  RRoverall0toTau <- 1 - VEoverall0toTau  
  
  # Use the following assumptions:
  # 1. The failure time T and the censoring time C are independent.
  # 2. T and C have exponential distributions in the hypothetical placebo arm
  #    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
  # 3. RRoverall = P(tau < T <= t|Z=1)/P(tau < T <= t|Z=0) for all t between
  #                tau and taumax (this will only approximately hold).
  
  # Calculate thetat and thetac on a monthly time-scale:
  thetat <- -log(1-risk0Mo0toTaumax)/taumax
  thetac <- -log(1-dropOutRate)/taumax
  
  # Calculate the number of subjects in the vaccine group at risk at month tau  
  # Logic: = Number enrolled * P(X > tau|Z=1)
  #        =    NvEnrolled   * [1-RRoverall0totau*P(T <= tau|Z=0)]*P(C > tau)
  #        =    NvEnrolled   * [1-RRoverall0totau*(1-exp(-thetat*tau))]*exp(-thetac*tau)
  N <- NvEnrolled*(1-RRoverall0toTau*pexp(tau,thetat))*(1-pexp(tau,thetac))
  # where pexp(t,theta) is the cumulative distribution function of an exponential random
  # variable with rate parameter theta.
  
  
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
  # nCases is notation in the paper. nCases is notation used in R program computepower
  nCases <- round(N*(term1+term2)/term3)
  
  # Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
  # Calculated the same as for term3 except using taumax instead of tau
  # nControls is notation in the paper. nControls is notation used in R program computepower
  nControls <- round(NvEnrolled*(1-RRoverall*pexp(taumax,thetat))*(1-pexp(taumax,thetac)))
  
  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax 
  # and with the biomarker measured (i.e., in Phase 2).
  # nCasesWithS is notation used in R program computepower
  nCasesWithS <- round(fracImmunog*nCases)
  
  return(list(N = N, nCases = nCases, nControls = nControls, nCasesWithS = nCasesWithS))
}


