#' Plot Power Curve against CoR Relative Risk for Trichotomous or Binary Biomarkers
#'
#' Plots the power to detect a trichotomous or binary CoR in vaccine recipients against different
#' values of the CoR effect size, RRt, which is the relative risk in the high biomarker
#' response subgroup versus the low biomarker response subgroup.
#'
#' @param outComputePower List or list of lists containing output from \code{\link{computePower}}, or character vector specifying the file(s) containing \code{\link{computePower}} output.
#' @param outDir Character vector specifying path(s) to output file(s), necessary if \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText Character vector specifying the entirety of the legend text. Order of the parameter values must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#'
#' @details If multiple levels are specified for the biomarker measurement error input parameters (i.e., for \code{Sens}/\code{Spec} or \code{rho}) in the
#' \code{\link{computePower}} function, only the first level is used to determine the \eqn{RR_t} values that are plotted on the x-axis.
#'
#' @return Plot displaying power vs. CoR relative risk
#'
#' @examples
#'
#' # Example scenario with trichotomous biomarker, where values of controlCaseRatio are varied
#'
#' # Set input parameters for computePower function
#' nCases <- 10
#' nControls <- 300
#' nCasesWithS <- 10
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' VElat0 <- seq(0, VEoverall, len=5)
#' VElat1 <- rep(VEoverall, 5)
#' Plat0 <- P0 <- 0.2
#' Plat2 <- P2 <- 0.6
#' sens <- spec <- 0.8
#' FP0 <- FN2 <- 0
#' M <- 50
#' alpha <- 0.05
#' biomType <- "trichotomous"
#'
#' # Output from computePower function is stored in an object as a list
#' pwr1 <- computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'                      controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                      Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2, VElat0=VElat0, VElat1=VElat1,
#'                      M=M, alpha=alpha, spec=spec, FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' controlCaseRatio <- 3
#' pwr2 <- computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'                      controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                      Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2, VElat0=VElat0, VElat1=VElat1,
#'                      M=M, alpha=alpha, spec=spec, FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' controlCaseRatio <- 1
#' pwr3 <- computePower(nCases=nCases, nControls=nControls, nCasesWithS=nCasesWithS,
#'                      controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                      Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2, VElat0=VElat0, VElat1=VElat1,
#'                      M=M, alpha=alpha, spec=spec, FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' # Set parameters for plotPowerTri function
#' # outComputePower is a list of lists containing outputs from the computePower function
#' outComputePower <- list(pwr1, pwr2, pwr3)
#' legendText <- paste0("controls:cases = ", c("5:1","3:1","1:1"))
#' plotPowerTri(outComputePower=outComputePower, legendText=legendText)
#'
#' \dontrun{
# # Output from computePower function is saved in RData files
# computePower(..., saveDir = "myDir", saveFile = "myFile.RData")

#' # outComputePower is a character vector specifying the files containing computePower output
#' # outDir is a character vector specifying the outComputePower file directories
#' outComputePower = c("myFile1", "myFile2", "myFile3")
#' outDir = rep("~/myDir", 3)
#' legendText <- paste0("controls:cases = ", c("5:1","3:1","1:1"))
#' plotPowerTri(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}} \code{\link{plotPowerCont}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotPowerTri <- function(outComputePower, outDir=NULL, legendText) {

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

  power <- pwr$power
  RRt <- pwr$RRt
  VElat0 <- pwr$VElat0
  VElat2 <- pwr$VElat2
  VElat1 <- rep(pwr$VEoverall, length(VElat0))
  alpha <- pwr$alpha
  Plat0 <- pwr$Plat0
  Plat2 <- pwr$Plat2
  biomType <- ifelse(Plat0 + Plat2 == 1, "Binary", "Trichotomous")

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
  power <- t(power)
  RRt <- t(RRt)

  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,oma=c(4,3,3,4))
  plot(RRt[,1],power[,1],ylim=c(0,1),type='n',xlab="",ylab="Power",axes=FALSE, cex.axis=1.2)
  box()

  m <- nrow(RRt)
  RRtgrid <- seq(RRt[1,1],RRt[m,1],len=5)
  axis(1,at=RRtgrid,labels=round(RRtgrid,2))
  axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0),labels=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0))
  tick1 <- which.min(abs(RRt[,1]-RRtgrid[1]))
  tick2 <- which.min(abs(RRt[,1]-RRtgrid[2]))
  tick3 <- which.min(abs(RRt[,1]-RRtgrid[3]))
  tick4 <- which.min(abs(RRt[,1]-RRtgrid[4]))
  tick5 <- which.min(abs(RRt[,1]-RRtgrid[5]))
  axis(1, at=RRtgrid, labels=round(c(VElat0[tick1], VElat0[tick2], VElat0[tick3], VElat0[tick4], VElat0[tick5]),2), line=1.5, tick=FALSE)
  if(biomType == "binary") {
    axis(1, at=RRtgrid, labels=round(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]),2), line=3, tick=FALSE)
  } else {
    axis(1, at=RRtgrid, labels=round(c(VElat1[tick1], VElat1[tick2], VElat1[tick3], VElat1[tick4], VElat1[tick5]),2), line=3, tick=FALSE)
    axis(1, at=RRtgrid, labels=round(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]),2), line=4.5, tick=FALSE)
  }
  mtext(expression(RR[t]), 1, line=1, at=-0.05, cex=1.2)
  mtext(expression(VE[0]^{lat}), 1, line=2.8, at=-0.05, cex=1.2)
  if(biomType == "binary") {
    mtext(expression(VE[2]^{lat}), 1, line=4.3, at=-0.05, cex=1.2)
  } else {
    mtext(expression(VE[1]^{lat}), 1, line=4.3, at=-0.05, cex=1.2)
    mtext(expression(VE[2]^{lat}), 1, line=5.8, at=-0.05, cex=1.2)
  }
  mtext(expression("Vaccine Group CoR Relative Risk"~ RR[t]),1, line=7, at=mean(RRtgrid), cex=1.2)

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRt[,1], power[,i], lty=i, col=colors[i], lwd=3)
  }

  abline(h=alpha/2,lty=3)

  legend(x="topright",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)

  title(bquote(paste("Power to Detect a "~.(biomType)~" CoR in Vaccine Recipients [2-sided ", alpha, "=", .(alpha), "]")))

}
