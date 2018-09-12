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
#' 
#' @details
#' This function performs simulations to calculate the power for testing whether a trichotomous or continuous biomarker
#' is a correlate of risk (CoR).
#'
#' @return N - the number of subjects in the vaccine group at risk at month tau
#'
#' @examples
#'
#' @import survival
#' @import osDesign
#' @export
computeN <- function(tau, taumax, VEoverall, NvEnrolled, risk0, dropOutRate) {
  
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
  
  return(N)
}


