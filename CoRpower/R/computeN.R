#' Sample Size Calculations for Correlates of Risk
#'
#' Calculates additional requisite design parameters pertaining to the target population of active treatment recipients observed to be
#' at risk at the biomarker sampling timepoint.
#'
#' @param Nrand the number of participants randomized to vaccine arm
#' @param tau the biomarker sampling timepoint
#' @param taumax the end of the follow-up time period
#' @param VEtauToTaumax the treatment (vaccine) efficacy level between \eqn{\tau} and \eqn{\tau_{max}}
#' @param VE0toTau the treatment (vaccine) efficacy between 0 and \eqn{\tau}
#' @param risk0 the placebo-group endpoint risk between \eqn{\tau} and \eqn{\tau_{max}}
#' @param dropoutRisk the risk of participant dropout between 0 and \eqn{\tau_{max}}
#' @param propCasesWithS the proportion of cases with measured biomarker S
#'
#' @details
#' The design parameters calculated by this function can be used as input values for the \code{\link{computePower}} function.
#' The calculations include options to account for participant dropout by specifying the \code{dropoutRisk} input parameter,
#' as well as for incomplete sample storage by specifying the \code{propCasesWithS} input parameter.
#'
#' The calculations use the following assumptions:
#' \enumerate{
#'   \item Failure time \eqn{T} and censoring time \eqn{C} are independent
#'   \item \eqn{T|Z=0} follows an exponential distribution with rate parameter \eqn{\theta_t} and \eqn{C|Z=0} follows an
#'      exponential distribution with rate parameter \eqn{\theta_c}
#'   \item RRtauToTaumax \eqn{ = P(T <= t|T> \tau, Z=1)/P(T <= t|T> \tau, Z=0)} for all \eqn{t} between tau and taumax (this will only approximately hold).
#' }
#' @return List with the following elements:
#' \itemize{
#'   \item N: the number of participants in the active treatment group that are at risk at \eqn{\tau}
#'   \item nCases: the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group
#'    \item nControls: the number of controls observed (or projected) to complete follow-up through \eqn{\tau_{max}} endpoint-free in the active treatment group
#'    \item nCasesWithS: the number of clinical endpoint cases observed (or projected) between \eqn{\tau} and \eqn{\tau_{max}} in the active treatment group with an available biomarker response
#' }
#'
#' @examples
#' Nrand = 4100
#' tau = 3.5
#' taumax = 24
#' VEtauToTaumax = 0.75
#' VE0toTau = 0.75/2
#' risk0 = 0.034
#' dropoutRisk = 0.1
#' propCasesWithS = 1
#' computeN(Nrand, tau, taumax, VEtauToTaumax, VE0toTau, risk0, dropoutRisk, propCasesWithS)
#'
#' @seealso \code{\link{computePower}}
#'
#' @importFrom stats pexp
#'
#' @export
computeN <- function(Nrand, tau, taumax, VEtauToTaumax, VE0toTau, risk0, dropoutRisk, propCasesWithS) {

  # With T the failure time, C the censoring time, X = min(T,C), and
  # Z vaccination status (vaccine or hypothetical placebo), we have
  # RRoverall = RRtauToTaumax = P(T <= taumax|T > tau, Z=1)/P(T <= taumax|T > tau, Z=0)
  RRoverall <- 1 - VEtauToTaumax

  # Hypothetical placebo-group endpoint risk between enrollment and time tau (default is in months):
  risk0Mo0toTaumax <- risk0*(1 + tau/taumax)

  # Specify relative risk in the early period between enrollment and time tau,
  # defined as RR0toTau = P(T <= tau|Z=1) / P(T <= tau|Z=0)
  RR0toTau <- 1 - VE0toTau

  # Use the following assumptions:
  # 1. The failure time T and the censoring time C are independent.
  # 2. T and C have exponential distributions in the hypothetical placebo arm
  #    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
  # 3. RRoverall = P(T <= t|T>tau, Z=1)/P(T <= t|T>tau, Z=0) for all t
  #    between tau and taumax (this will only approximately hold)

  # Calculate thetat and thetac on the same time-scale as that of tau and taumax:
  thetat <- -log(1-risk0Mo0toTaumax)/taumax
  thetac <- -log(1-dropoutRisk)/taumax

  # Calculate the number of subjects in the vaccine group at risk at tau
  # Logic: = Number enrolled * P(X > tau|Z=1)
  #        =    Nrand   * [1-RR0toTau*P(T <= tau|Z=0)]*P(C > tau)
  #        =    Nrand   * [1-RR0toTau*(1-exp(-thetat*tau))]*exp(-thetac*tau)
  N <- Nrand*(1-RR0toTau*pexp(tau,thetat))*(1-pexp(tau,thetac))
  # where pexp(t,theta) is the cumulative distribution function of an exponential random
  # variable with rate parameter thetat.


  # Number of observed rotavirus endpoints between tau and taumax in vaccinees:
  # Logic: = Number vaccinees at-risk at tau * P(T<=taumax,T<=C|X>tau,Z=1)
  #        = N*\int_0^infty P(T<=min(c,taumax),C=c|X>tau,Z=1)dc
  #        = N*\int_tau^infty P(tau<T<=min(c,taumax),C=c|Z=1)dc / P(X>tau|Z=1)
  #        = N*{\int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc}/P(X>tau|Z=1).

  # Now, the above integral in { } can be written as
  #   \int_tau^infty P(tau<T<=min(c,taumax)|Z=1)P(C=c)dc
  # = \int_tau^taumax P(tau<T<=c|Z=1)P(C=c)dc + P(tau<T<=taumax|Z=1)\int_taumax^infty P(C=c)dc
  # = \int_tau^taumax {P(T<=c|T>tau,Z=1)P(T>tau|Z=1)P(C=c)} + P(T<=taumax|T>tau,Z=1)\int_taumax^infty P(C=c)dc
  # = RRoverall*(1-RR0toTau*pexp(tau, thetat))*\int_tau^taumax{P(tau<T<=c|Z=0)/P(T>tau|Z=0)*thetac*exp(-thetac*c)dc} (term1)
  # + RRoverall*(1-RR0toTau*pexp(tau, thetat))*P(tau<T<=c|Z=0)/P(T>tau|Z=0)*[1 - pexp(taumax, thetac)] (term2)

  # Term (1) is calculated to be:
  #   RRoverall*(1-RR0toTau*pexp(tau, thetat))/P(T>tau|Z=0)*[exp(-thetat*tau)*(pexp(taumax,thetac)-pexp(tau,thetac))
  # - thetac*\int_tau^taumax exp(-(thetat+thetac)c)dc]
  # = RRoverall*(1-RR0toTau*pexp(tau, thetat))*[(pexp(taumax,thetac)-pexp(tau,thetac))
  # - thetac/((thetac+thetat)*exp(-thetat*tau))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))]

  # Term (2) is calculated to be:
  #  RRoverall*(1-RR0toTau*pexp(tau, thetat))*[(pexp(taumax, thetat)-pexp(tau, thetat))/exp(-thetat*tau)]*[1 - pexp(taumax, thetac)]

  # Lastly, the denominator (call it term3) equals:
  # P(X>tau|Z=1) = [1-RR0toTau*P(T<=tau|Z=0)]P(C>tau)
  #              = [1-RR0toTau*pexp(tau,thetat)][1-pexp(tau,thetac)].

  # Putting it together
  term1 <- RRoverall*(1-RR0toTau*pexp(tau, thetat))*((pexp(taumax,thetac)-pexp(tau,thetac))
                      - (thetac/((thetac+thetat)*exp(-thetat*tau)))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))

  term2 <- RRoverall*(1-RR0toTau*pexp(tau, thetat))*((pexp(taumax,thetat)-pexp(tau,thetat))/exp(-thetat*tau))*(1-pexp(taumax,thetac))
  term3 <- (1-RR0toTau*pexp(tau,thetat))*(1-pexp(tau,thetac))

  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax.
  # nCases is notation used in the paper and in the R program computepower
  nCases <- round(N*(term1+term2)/term3)

  # Number of subjects in the vaccine group at-risk at tau without the clinical event (controls) by taumax.
  # Calculated the same as for term3 except using taumax instead of tau.
  # nControls is notation used in the paper and in the R program computepower
  nControls <- round(Nrand*(1-RRoverall*pexp(taumax,thetat))*(1-pexp(taumax,thetac)))

  # Number of subjects in the vaccine group at-risk at tau with the clinical event (cases) by taumax and with the biomarker measured.
  # nCasesWithS is notation used in R program computepower
  nCasesWithS <- round(propCasesWithS*nCases)

  return(list(N = round(N), nCases = nCases, nControls = nControls, nCasesWithS = nCasesWithS))
}


