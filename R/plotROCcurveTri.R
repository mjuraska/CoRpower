# Internal function

computesensspecFPFN <- function(sigma2obs,rho,Plat0,Plat2,P0,P2) {
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
  #   P0: Probability of low biomarker response. May be scalar or vector specifying grid of probabilities.
  #       If vector, should include one value equal to Plat0 and values straddling either side.
  #   P2: Probability of high biomarker response. May be scalar or vector specifying grid of probabilities.
  #       If vector, should include one value equal to Plat2 and values straddling either side.
  #
  # Returns:
  #   Matrix with each row corresponding to one set of values for Plat0, Plat2, P0, and P2,
  #   and with the following columns: thetaloVE*, thetahiVE*, Plat0, Plat1, Plat2, P0, P2,
  #   taulosolution**, tauhisolution**, sens, spec, FP0, FP1, FN2, and FN1.
  #   If dichotomous biomarker, returns 0's for FP0, FP1, FN2, and FN1, which are irrelavant.
  #
  #   *Let X* denote the true latent subgroups. X* = 2 if X* > thetahiVE, X* = 0 if X <= thetaloVE, and
  #   X* = 1 if X* is in between thetaloVE and thetahiVE, for fixed thetaloVE and thetahiVE that are solved for.
  #   **Let S* denote a trichotomous biomarker. S* = 2 if S* > tauhisolution and S* = 0 if S* <= taulosolution,
  #   and S* = 1 if S* is in between taulosolution and tauhisolution, for fixed taulosolution and tauhisolution
  #   that are solved for.

  sigma2e <- (1-rho)*sigma2obs
  sigma2tr <- rho*sigma2obs
  thetahiVE <- qnorm(1-Plat2)*sqrt(sigma2tr)
  thetaloVE <- qnorm(Plat0)*sqrt(sigma2tr)

  Plat1 <- 1 - Plat0 - Plat2
  m <- length(P2)
  Sens <- rep(1,m)
  Spec <- rep(1,m)
  FP1 <- rep(0,m)
  FP2 <- rep(0,m)
  FN1 <- rep(0,m)
  FN2 <- rep(0,m)
  tauhisolution <- rep(0,m)
  taulosolution <- rep(0,m)

  if (rho < 1) {
    # Stochastic integration
    X <- rnorm(10000,0,sqrt(sigma2tr))
    S <- X + rnorm(10000,0,sqrt(sigma2e))

    Phi <- sum(X>thetahiVE)/length(X)
    Plo <- sum(X<=thetaloVE)/length(X)
    Pmed <- 1 - Phi - Plo

    for (l in 1:m) {
      # Find the cut points tauhi and taulo by solving the following equations:
      #   0 = sensvec*Plat2 + FP1vec*Plat1 + FP0vec*Plat0 - P2  (f2 below; eqn 8 in manuscript)
      #   0 = specvec*Plat0 + FN1vec*Plat1 + FN2vec*Plat2 - P0  (f0 below; eqn 7 in manuscript)
      # where
      #   sensvec <- (sum(S>tauhi & X > thetahiVE)/length(S))/Phi
      #   specvec <- (sum(S<=taulo & X <= thetaloVE)/length(S))/Plo
      # if dichotomous biomarker,
      #   FP1vec <- 0
      #   FP0vec <- 0
      #   FN2vec <- 0
      #   FN1vec <- 0
      # if trichotomous biomarker,
      #   FP1vec <- (sum(S>tauhi & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed
      #   FP0vec <- (sum(S>tauhi & X <= thetaloVE)/length(S))/Plo
      #   FN2vec <- (sum(S<=taulo & X > thetahiVE)/length(S))/Phi
      #   FN1vec <- (sum(S<=taulo & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed

      if (Pmed==0){  # dichotomous
        f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE)/length(S))/Phi)*Plat2 - P2[l]
        f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE)/length(S))/Plo)*Plat0 - P0[l]
      } else {  # trichotomous
        f2 <- function(tauhi) ((sum(S>tauhi & X > thetahiVE)/length(S))/Phi)*Plat2 +
          ((sum(S>tauhi & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed)*Plat1 +
          ((sum(S>tauhi & X <= thetaloVE)/length(S))/Plo)*Plat0 - P2[l]
        f0 <- function(taulo) ((sum(S<=taulo & X <= thetaloVE)/length(S))/Plo)*Plat0 +
          ((sum(S<=taulo & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed)*Plat1 +
          ((sum(S<=taulo & X > thetahiVE)/length(S))/Phi)*Plat2 - P0[l]
      }

      tauhisol <- uniroot(f2, interval=c(-2.5, 2.5))$root
      taulosol <- uniroot(f0, interval=c(-2.5, 2.5))$root
      tauhisolution[l] <- tauhisol
      taulosolution[l] <- taulosol

      Sens[l] <- sum(S>tauhisolution[l] & X > thetahiVE)/sum(X>thetahiVE)
      Spec[l] <- sum(S<=taulosolution[l] & X <= thetaloVE)/sum(X<=thetaloVE)
      if (Pmed==0) {
        FP2[l] <- 0
        FP1[l] <- 0
        FN2[l] <- 0
        FN1[l] <- 0 }
      if (Pmed > 0) {
        FP2[l] <- (sum(S>tauhisolution[l] & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed
        FP1[l] <- (sum(S>tauhisolution[l] & X <= thetaloVE)/length(S))/Plo
        FN1[l] <- (sum(S<=taulosolution[l] & X > thetahiVE)/length(S))/Phi
        FN2[l] <- (sum(S<=taulosolution[l] & X > thetaloVE & X <= thetahiVE)/length(S))/Pmed
      }
    }
  }

  ans <- cbind(rep(thetaloVE,m),rep(thetahiVE,m),rep(Plat0,m),rep(Plat1,m),rep(Plat2,m),P0,P2,
               taulosolution,tauhisolution,Sens,Spec,FP1,FP2,FN1,FN2)
  return(ans)
}

#' Plotting of ROC Curves for Trichotomous Biomarkers
#'
#' Plots the receiver operating characteristic (ROC) curve displaying sensitivity and specificity for a range of \code{P2} and \code{P0} values,
#' four values of \code{rho}, and four values of \code{Plat2}. Illustrates how different levels of measurement error \code{rho} map to sensitivity
#' and specificity, depending on the value of \code{Plat2}. This funciton is used to create Figure 1 in the Supplementary Material of
#' [Gilbert, Janes, and Huang (2016). “Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials.”]
#'
#' @param Plat0 a numeric value specifying the prevalence of the latent lower protected subgroup for a dichotomous or trichotomous biomarker
#' @param Plat2 a numeric vector of length four specifying the prevalences of the latent higher protected subgroup for a dichotomous or trichotomous biomarker
#' @param P0 a numeric vector specifying a grid of probabilities of low biomarker response for a dichotomous or trichotomous biomarker.
#' @param P2 a numeric vector specifying a grid of probabilities of high biomarker response for a dichotomous or trichotomous biomarker.
#' @param rho a numeric vector of length four specifying distinct protection-relevant fractions of \code{sigma2obs}.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' Plat0 <- 0.2
#' Plat2 <- c(0.2, 0.3, 0.4, 0.5)
#' P0 <- seq(0.90, 0.10, len=25)
#' P2 <- seq(0.10, 0.90, len=25)
#' rho <- c(1, 0.9, 0.7, 0.5)
#' plotROCcurveTri(Plat0 = Plat0, Plat2 = Plat2, P0 = P0, P2 = P2, rho = rho)
#'
#' @export
plotROCcurveTri <- function(Plat0, Plat2, P0, P2, rho) {
  if(length(rho)!= 4 | any(rho > 1) | any(rho < 0)) {
    stop("rho must be a numeric vector of 4 elements, with each element's value between 0 and 1")
  } else if(length(Plat2)!= 4 | any(Plat2 > 1) | any(Plat2 < 0)) {
    stop("Plat2 must be a numeric vector of 4 elements, with each element's value between 0 and 1")
  }

  par(cex.axis=1,cex.lab=1,cex.main=1.2,mar=c(4.1,4.1,4.1,4.1),mfrow=c(2,2),oma=c(3,3,5,5),las=1)
  for(i in 1:length(Plat2)) {
    Plat2VE <- Plat2[i]
    ss1 <- computesensspecFPFN(1, rho[1],Plat0,Plat2VE,P0,P2)
    ss2 <- computesensspecFPFN(1, rho[2],Plat0,Plat2VE,P0,P2)
    ss3 <- computesensspecFPFN(1, rho[3],Plat0,Plat2VE,P0,P2)
    ss4 <- computesensspecFPFN(1, rho[4],Plat0,Plat2VE,P0,P2)

    plot(1-ss1[,11],ss1[,10],xlim=c(0,0.5),ylim=c(0.5,1),type='n',axes=FALSE,xlab="1-Specificity = 1-P(S=0|X=0)",ylab="Sensitivity = P(S=2|X=2)")
    axis(1)
    axis(2)
    lines(1-ss1[,11],ss1[,10],lty=1,col="blue",lwd=4)
    lines(1-ss2[,11],ss2[,10],lty=2,col="orange",lwd=4)
    lines(1-ss3[,11],ss3[,10],lty=3,col="red",lwd=4)
    lines(1-ss4[,11],ss4[,10],lty=4,col="black",lwd=4)
    legend(x="bottomright",legend=c(as.expression(bquote(rho==.(rho[1]))),as.expression(bquote(rho==.(rho[2]))),
                                    as.expression(bquote(rho==.(rho[3]))),as.expression(bquote(rho==.(rho[4])))),
           lty=c(1,2,3,4),col=c("blue","orange","red","black"),lwd=3)
    title(bquote(P[2]^{lat}==.(Plat2VE)), cex=1.1)
  }

  mtext(bquote("ROC Curve of a Trichotomous Marker: "~.(P2[1]*100)~"% - "~.(P2[length(P2)]*100)~"%  ("~.(P0[1]*100)~"% - "~.(P0[length(P0)]*100)~"%) Vaccinees with S=2 (S=0)"),outer=T,cex=1.2)
}



