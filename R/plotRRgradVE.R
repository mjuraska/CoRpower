#' Plotting of the Ratio of Relative Risks for Higher/Lower Latent Subgroups against Correlate of Risk Effect Size for Trichotomous Biomarkers
#'
#' Plots the ratio of relative risks for the higher and lower latent subgroups (on the y-axis) versus the correlate of risk effect size (on the x-axis) in the active treatment group for a trichotomous biomarker.
#' The correlate of risk effect size is quantified as the relative risk ratio of the clinical endpoint comparing subgroups of active treatment recipients with high and low biomarker response.
#'
#' @param outComputePower either a list or list of lists containing output from \code{\link{computePower}} or a character vector specifying the \code{.RData} file(s) containing \code{\link{computePower}} output
#' @param outDir a character vector specifying path(s) to output \code{.RData} file(s), necessary if \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText a character vector specifying the entirety of the legend text. The order of the elements (i.e., parameter values) must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#'
#' @details
#' When \code{rho} is varied, this plot shows how the relationship between the correlate of risk effect size and the relative risks for the higher and lower latent subgroups
#' changes for different values of \code{rho}. The ratio of relative risks for the higher and lower latent subgroups is a relative vaccine efficacy parameter. When \code{rho=1},
#' a correlate of risk in the vaccine group is equivalent to the relative vaccine efficacy parameter, whereas for imperfectly measured biomarkers with \code{rho<1},
#' the correlate of risk effect size is closer to the null than the relative vaccine efficacy parameter is.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' # Example scenario with trichotomous biomarker, where values of rho are varied
#'
#' # Set input parameters for computePower function
#' nCases <- 10
#' nControls <- 300
#' nCasesWithS <- 10
#' controlCaseRatio <- 3
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' VElat0 <- seq(0, VEoverall, len=10)
#' VElat1 <- rep(VEoverall, 10)
#' Plat0 <- P0 <- 0.2
#' Plat2 <- P2 <- 0.6
#' M <- 20
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.7, 0.4)
#' biomType <- "trichotomous"
#'
#' # Output from computePower function is stored in an object as a list
#' pwr <- computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'                     controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                     Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2, VElat0=VElat0, VElat1=VElat1,
#'                     M=M, alpha=alpha, sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#' # Set parameters for plotPowerCont function
#' # outComputePower is a list containing output from the computePower function
#' outComputePower <- pwr
#' legendText <- paste0("rho = ", c(1, 0.7, 0.4))
#' plotRRgradVE(outComputePower=outComputePower, legendText=legendText)
#'
#' \dontrun{
#' # Output from computePower function is saved in an RData file
#' computePower(..., saveDir = "myDir", saveFile = "myFile.RData")

#' # outComputePower is a character string specifying the file containing the computePower output
#' # outDir is a character string specifying the outComputePower file directory
#' outComputePower = "myFile"
#' outDir = "~/myDir"
#' legendText <- paste0("rho = ", c(1, 0.7, 0.4))
#' plotRRgradVE(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}} \code{\link{plotPowerTri}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotRRgradVE <- function(outComputePower, outDir=NULL, legendText) {

  multiple <- TRUE
  if(any(sapply(outComputePower, is.list))) {  # check if list of lists
    pwr <- outComputePower[[1]]  # load first output list
  } else if(is.list(outComputePower)) {  # check if single list
    pwr <- outComputePower
    multiple <- FALSE
  } else if(is.character(outComputePower) & is.null(outDir)) {  # check outDir is specified
    stop("outComputePower is a character vector so outDir needs to be specified")
  } else if(is.character(outComputePower)) {  # check if character
    load(paste0(file.path(outDir[1], outComputePower[1]),".RData"))  # load first output list
  } else {
    stop("outComputePower must be of type list or character")
  }

  power <- t(pwr$power)
  RRt <- t(pwr$RRt)
  alpha <- pwr$alpha
  VElat2 <- pwr$VElat2
  VElat0 <- pwr$VElat0
  RRlat2 <- 1 - VElat2
  RRlat0 <- 1 - VElat0
  ratio <- RRlat2/RRlat0

  if(multiple==TRUE) {
    for(i in 2:length(outComputePower)) {
      if(is.list(outComputePower)) {
        pwr <- outComputePower[[i]]
      } else {
        load(paste0(file.path(outDir[i], outComputePower[i]),".RData"))
      }
      addPower <- pwr$power
      power <- rbind(power, addPower)
    }
  }

  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,oma=c(4,3,4,4), mar=c(4,5,2,2))
  plot(RRt[,1],ratio,xlim=c(0,1),ylim=c(0,1),type='n',axes=FALSE,ylab=expression(RR[2]^{lat} / RR[0]^{lat}),xlab=expression("Vaccine Group CoR Relative Risk"~ RR[t]==risk[1](2)/risk[1](0)))
  axis(1)
  axis(2)
  box()

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRt[,i], ratio, lty=i, col=colors[i], lwd=3)
  }
  abline(0,1,lty=3)

  title(expression("RR Ratio in the Higher/Lower Latent Subgroups vs. CoR Relative Risk"~RR[t]))

  legend(x="topleft",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)

}
