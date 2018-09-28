#' Plot Vaccine Efficacy Curves for Different CoR Relative Risks for Continuous Biomarkers
#'
#' Plots the vaccine efficacy (VE) curve for the true biomarker X*=x* for eight different values of the true CoR relative risk, \eqn{RR_c (\rho=1)},
#' in vaccine recipients and the lowest vaccine efficacy level for the true biomarker, \eqn{VE_lowest}.
#' All curves assume \eqn{\rho=1}, and VE ranges from 0 to 1. The legend is completely determined by the function.
#'
#' @param outComputePower List containing output from \code{\link{computePower}}, or character string specifying the file containing \code{\link{computePower}} output. Must be one list or character string.
#' @param outDir Character string specifying path to output file, necessary if \code{outComputePower} is a character string. Default is \code{NULL}. Must be one character string; cannot take a character vector.
#'
#' @details
#' \code{\link{computePower}} function input parameter \code{VElowest} must have length >= 8
#' for all eight scenarios to have unique RRc and VElowest. Otherwise, only \code{length(VElowest)} unique
#' VE curves will be displayed.
#'
#' When interpreting the output of the function, the null hypothesis \eqn{H_0: CoR RR_c=1} (expression 16 in the manuscript)
#' corresponds to a flat curve where VE(x*) = VE for all x*. Increasing departures from the null hypothesis correspond
#' to increasingly variable and steep VE curves. The output assumes \eqn{risk_0} is constant for all x* and s* and there is
#' no measurement error (\eqn{\rho=1}), and one can see that when this is the case, an association of the biomarker with
#' infection risk in the vaccine group (a CoR) is equivalent to an association of the biomarker with VE.
#'
#' The function's plot can also be interpreted in conjunction with the output of the \code{\link{plotPowerCont}} function by
#' matching the CoR relative risk in the two plots and examining power compared to VE. This sheds light on the importance
#' of overall VE on power and further enables correlates of risk results to be interpreted in terms of
#' potential correlates of efficacy/protection.
#'
#' @return Plot displaying VE curves for various CoR relative risk and lowest VE scenarios
#'
#' @examples
#'
#' # Example scenario with continuous biomarker, where values of rho are varied
#'
#' # Set input parameters for computePower function
#' nCases <- 10
#' nControls <- 300
#' nCasesWithS <- 10
#' controlCaseRatio <- 3
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' PlatVElowest <- 0.2
#' VElowest <- seq(0, VEoverall, len=8)
#' M <- 13
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.7, 0.4)
#' biomType <- "continuous"
#'
#' # Output from computePower function is stored in an object as a list
#' pwr <- computePower(nCases=nCases, nCasesWithS=nCasesWithS, nControls=nControls,
#'                     controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                     PlatVElowest=PlatVElowest, VElowest=VElowest, M=M, alpha=alpha,
#'                     sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#' # Set parameters for plotPowerCont function
#' # outComputePower is a list containing output from the computePower function
#' outComputePower <- pwr
#' plotVElatCont(outComputePower=outComputePower)
#'
#' \dontrun{
#' # Output from computePower function is saved in an RData file
#' computePower(..., saveDir = "myDir", saveFile = "myFile.RData")

#' # outComputePower is a character string specifying the file containing the computePower output
#' # outDir is a character string specifying the outComputePower file directory
#' outComputePower = "myFile"
#' outDir = "~/myDir"
#' plotVElatCont(outComputePower, outDir=outDir)
#' }
#'
#' @seealso \code{\link{computePower}} \code{\link{plotPowerCont}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotVElatCont <- function(outComputePower, outDir=NULL) {
  if(any(sapply(outComputePower, is.list)) | length(outDir)>1) {
    stop("outComputePower must be a single list, not a list of lists, and outDir must be of length 1")
  } else if(is.list(outComputePower)) {
    pwr <- outComputePower
  } else if(is.character(outComputePower) & is.null(outDir)) {
    stop("outComputePower is a character vector so outDir needs to be specified")
  } else if(is.character(outComputePower)) {
    load(paste0(file.path(outDir[1], outComputePower[1]),".RData"))
  } else {
    stop("outComputePower must be of type list or character")
  }

  RRc <- pwr$RRc
  alpha <- pwr$alpha
  sigma2obs <- pwr$sigma2obs
  PlatVElowest <- pwr$PlatVElowest
  VElowest <- pwr$VElowest
  VEoverall <- pwr$VEoverall
  risk0 <- pwr$risk0

  rho <- 1
  nu <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  o <- length(VElowest)

  betaLat <- rep(NA,o)
  alphaLat <- rep(NA,o)
  for (l in 1:o) {
    # find solutions alphalat and betalat by solving eqn (4) in Appendix B
    risk1latnu <- (1-VElowest[l])*risk0

    alphaLat[l] <- uniroot(alphaLatEqn, lower=-10, upper=10, nu=nu, risk1latnu=risk1latnu, sigma2obs=sigma2obs, VEoverall=VEoverall, PlatVElowest=PlatVElowest, risk0=risk0)$root

    # Second solve for betalat:
    D <- risk1latnu
    betaLat[l] <- (log(D/(1-D)) - alphaLat[l])/nu[1]
  }

  svect <- seq(-3,3,len=300)

  # Compute VE(s_1) vs s_1 for each fixed RRc value:
  inds <- round(seq(1,o,len=8))
  fixedBetaLat <- betaLat[inds]
  fixedRRc <- RRc[inds]

  VEcurverr <- list()
  for(i in 1:8) {
    linpart <- alphaLat[inds[9-i]] + fixedBetaLat[9-i]*svect
    VEcurverrTemp <- 1-(exp(linpart)/(1+exp(linpart)))/risk0
    VEcurverrTemp[svect <= nu] <- VElowest[inds[9-i]]
    VEcurverr[[i]] <- VEcurverrTemp
  }

  ylims <- c(-19,0.95)
  newylims <- log10(1-(ylims))


  par(cex.axis=1.4,cex.lab=1.4,cex.main=1.35,mar=c(4,5,4,4),oma=c(3,0,3,0),las=1)
  plot(svect,VEcurverr[[1]],xlim=c(-3,3),ylim=c(-0.05,1),type='n',xlab=expression("True Biomarker X*=x* with "~rho==1~"in Vaccinees"), ylab=expression(VE[paste(x,"*")]^{lat}), axes=FALSE)
  axis(1)
  axis(2)

  colors <- c("brown","blue","dark grey","orange","red","purple","green","black")
  for(i in 1:8) {
    lines(svect, VEcurverr[[i]], lty=i, col=colors[i], lwd=4)
  }

  legend(x="bottomright",legend=c(as.expression(bquote("CoR RR=1.00,   "~VE[lowest]==.(round(VElowest[inds[8]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[7],2))~", "~VE[lowest]==.(round(VElowest[inds[7]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[6],2))~", "~VE[lowest]==.(round(VElowest[inds[6]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[5],2))~", "~VE[lowest]==.(round(VElowest[inds[5]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[4],2))~", "~VE[lowest]==.(round(VElowest[inds[4]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[3],2))~", "~VE[lowest]==.(round(VElowest[inds[3]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[2],2))~", "~VE[lowest]==.(round(VElowest[inds[2]],2)))),
                                  as.expression(bquote("CoR RR="~.(round(fixedRRc[1],2))~",      "~VE[lowest]==.(round(VElowest[inds[1]],2))))),
         lty=c(1,2,3,4,5,6,7,8),col=c("brown","blue","dark grey","orange","red","purple","green","black"),lwd=2, cex=1.2)
  title("VE Curves for Different CoR Relative Risks per 1 SD in Vaccinees")
  text(2, 0.46, bquote(atop(VE[overall]==.(round(VEoverall,2))~", "~risk[0]==.(round(risk0,3)),
                            ~P[lowestVE]^{lat}==.(round(PlatVElowest,2)))),cex=1.3)

}

