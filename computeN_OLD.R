#' Sample size calculations for assessing biomarkers as correlates of risk (CoRs) 
#'
#' Calculates projected sample sizes at design stage of trial assessing biomarkers as correlates of risk (CoRs) [Gilbert, Janes, and Huang (2015).
#' 
#' @param Nrand Number of participants randomized to vaccine arm
#' @param tau Biomarker sampling timepoint
#' @param taumax End of follow-up time period
#' @param VEtauToTaumax VE between 'tau' and 'taumax'
#' @param VE0toTau VE between 0 and 'tau'
#' @param risk0 Placebo-group endpoint risk between 'tau' and 'taumax'
#' @param dropoutRisk Dropout risk between 0 and 'taumax'
#' @param propCasesWithS Proportion of cases with measured S
#' 
#' @details
#' This function calculates projected sample sizes and accounts for dropout and incomplete sample storage. 
#'
#' @return List with the following elements: 
#'   N: the number of subjects in the vaccine group at risk at month tau
#'   nCases: Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax.
#    nControls: Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
#    nCasesWithS: Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax and with the biomarker measured
#'     
#'
#' @examples
#' Nrand = 4100,          
#' tau = 3.5,             
#' taumax = 24,           
#' VEtauToTaumax = 0.75, 
#' VE0toTau = 0.75/2,    
#' risk0 = 0.034,         
#' dropoutRisk = 0.1,     
#' propCasesWithS = 1    
#' computeN(Nrand, tau, taumax, VEtauToTaumax, VE0toTau, risk0, dropoutRisk, propCasesWithS)
#' 
#'
#' @import survival
#' @import osDesign
#' @export
computeN <- function(Nrand, tau, taumax, VEtauToTaumax, VE0toTau, risk0, dropoutRisk, propCasesWithS) {
  
  # With T the failure time, C the censoring time, X = min(T,C), and 
  # Z vaccination status (vaccine or hypothetical placebo), we have
  # RRoverall = P(tau < T <= taumax|Z=1)/P(tau < T <= taumax|Z=0) 
  RRoverall <- 1 - VEtauToTaumax
  
  # Cumulative failure rate in the hypothetical placebo arm between Month 0 and 24:
  risk0Mo0toTaumax <- risk0*(1 + tau/taumax)
  
  # Specify VE in the early period between enrollment and time tau 
  RRoverall0toTau <- 1 - VE0toTau
  
  # Use the following assumptions:
  # 1. The failure time T and the censoring time C are independent.
  # 2. T and C have exponential distributions in the hypothetical placebo arm
  #    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
  # 3. RRoverall = P(tau < T <= t|Z=1)/P(tau < T <= t|Z=0) for all t between
  #                tau and taumax (this will only approximately hold).
  
  # Calculate thetat and thetac on a monthly time-scale:
  thetat <- -log(1-risk0Mo0toTaumax)/taumax
  thetac <- -log(1-dropoutRisk)/taumax
  
  # Calculate the number of subjects in the vaccine group at risk at month tau  
  # Logic: = Number enrolled * P(X > tau|Z=1)
  #        =    Nrand   * [1-RRoverall0totau*P(T <= tau|Z=0)]*P(C > tau)
  #        =    Nrand   * [1-RRoverall0totau*(1-exp(-thetat*tau))]*exp(-thetac*tau)
  N <- Nrand*(1-RRoverall0toTau*pexp(tau,thetat))*(1-pexp(tau,thetac))
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
  nControls <- round(Nrand*(1-RRoverall*pexp(taumax,thetat))*(1-pexp(taumax,thetac)))
  
  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax 
  # and with the biomarker measured (i.e., in Phase 2).
  # nCasesWithS is notation used in R program computepower
  nCasesWithS <- round(propCasesWithS*nCases)
  
  return(list(N = round(N), nCases = nCases, nControls = nControls, nCasesWithS = nCasesWithS))
}


