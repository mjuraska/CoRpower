#' Plotting of Power versus Correlate of Risk Effect Size for Dichotomous and Trichotomous Biomarkers
#'
#' Plots power (on the y-axis) to detect a correlate of risk effect size (on the x-axis) in the active treatment group for a dichotomous or trichotomous biomarker. The correlate of risk effect size is quantified as
#' the relative risk ratio of the clinical endpoint comparing subgroups of active treatment recipients with high and low biomarker response.
#'
#' @param outComputePower either a list of lists containing output from \code{\link{computePower}} or a character vector specifying the \code{.RData} file(s) containing \code{\link{computePower}} output
#' @param outDir a character vector specifying path(s) to output \code{.RData} file(s), necessary if \cr
#' \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText a character vector specifying the entirety of the legend text. The order of the elements (i.e., parameter values) must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#'
#' @details If multiple levels are specified for the biomarker measurement error input parameters (i.e., for \code{sens}/\code{spec} or \code{rho}) in \code{\link{computePower}}, only the first level is used to determine
#' the \eqn{RR_t} values that are plotted on the x-axis.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' # Example scenario with trichotomous biomarker, where values of controlCaseRatio are varied
#'
#' # Set input parameters for computePower function
#' nCasesTx <- 10
#' nControlsTx <- 300
#' nCasesTxWithS <- 10
#' controlCaseRatio <- c(5,3)
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
#' # Output from computePower function is stored in an object as a list of lists
#' pwr <- computePower(nCasesTx=nCasesTx, nControlsTx=nControlsTx, nCasesTxWithS=nCasesTxWithS,
#'                      controlCaseRatio=controlCaseRatio, risk0=risk0,
#'                      VEoverall=VEoverall, Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2,
#'                      VElat0=VElat0, VElat1=VElat1, M=M, alpha=alpha, spec=spec,
#'                      FP0=FP0, sens=sens, FN2=FN2, biomType=biomType)
#'
#' # Set parameters for plotPowerTri function
#' # outComputePower is a list of lists containing outputs from the computePower function
#' outComputePower <- pwr
#' legendText <- paste0("controls:cases = ", c("5:1","3:1"))
#' plotPowerTri(outComputePower=outComputePower, legendText=legendText)
#'
#' \dontrun{
# # Output from computePower function is saved in RData files
# computePower(..., saveDir = "myDir", saveFile = "myFile.RData")

#' # outComputePower is a character vector specifying the files containing computePower output
#' # outDir is a character vector specifying the outComputePower file directories
#' outComputePower = c("myFile_controlCaseRatio_5.RData", "myFile_controlCaseRatio_3.RData")
#' outDir = rep("~/myDir", 2)
#' legendText <- paste0("controls:cases = ", c("5:1","3:1"))
#' plotPowerTri(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}}, \code{\link{plotPowerCont}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotPowerTri <- function(outComputePower, outDir=NULL, legendText) {

  if(any(sapply(outComputePower, is.list))) {  # check if list of lists
    pwr <- outComputePower[[1]]  # load first output list
  } else if(is.character(outComputePower) & is.null(outDir)) {  # check outDir is specified
    stop("outComputePower is a character vector so outDir needs to be specified")
  } else if(is.character(outComputePower)) {  # check if character
    fileName <- outComputePower[[1]]
    if (substr(fileName, start = nchar(fileName) - 5, stop = nchar(fileName)) != ".RData") {
      stop("File name(s) in outComputePower must include .RData")
    }
    load(file.path(outDir[1], fileName))  # load first output list
  } else {
    stop("outComputePower must be of type list or character")
  }

  power <- pwr$power
  RRt <- pwr$RRt
  VEoverall <- pwr$VEoverall
  VElat0 <- pwr$VElat0
  VElat2 <- pwr$VElat2
  alpha <- pwr$alpha
  Plat0 <- pwr$Plat0
  Plat2 <- pwr$Plat2
  Plat1 <- 1 - Plat0 - Plat2
  VElat1 <- (VEoverall - Plat0*VElat0 - Plat2*VElat2)/Plat1
  biomType <- ifelse(Plat1 == 0, "Dichotomous", "Trichotomous")

  if(length(outComputePower) > 1) {
    for(i in 2:length(outComputePower)) {
      if(is.list(outComputePower)) {
        pwr <- outComputePower[[i]]
      } else {
        load(paste0(file.path(outDir[i], outComputePower[i]),".RData"))
      }
      addPower <- pwr$power
      power <- cbind(power, addPower)
    }
  }
  
  power <- as.matrix(power)

  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,oma=c(4,3,3,4))
  plot(RRt,power[,1],ylim=c(0,1),type='n',xlab="",ylab="Power",axes=FALSE, cex.axis=1.2)
  box()

  RRtgrid <- seq(RRt[1],RRt[length(RRt)],len=5)
  axis(1, at=RRtgrid, labels=format(RRtgrid, digits=2, nsmall=2))
  axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0),labels=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0))
  tick1 <- which.min(abs(RRt-RRtgrid[1]))
  tick2 <- which.min(abs(RRt-RRtgrid[2]))
  tick3 <- which.min(abs(RRt-RRtgrid[3]))
  tick4 <- which.min(abs(RRt-RRtgrid[4]))
  tick5 <- which.min(abs(RRt-RRtgrid[5]))
  axis(1, at=RRtgrid, labels=format(c(VElat0[tick1], VElat0[tick2], VElat0[tick3], VElat0[tick4], VElat0[tick5]), digits=2, nsmall=2), line=1.5, tick=FALSE)
  if(biomType == "Dichotomous") {
    axis(1, at=RRtgrid, labels=format(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]), digits=2, nsmall=2), line=3, tick=FALSE)
  } else {
    axis(1, at=RRtgrid, labels=format(c(VElat1[tick1], VElat1[tick2], VElat1[tick3], VElat1[tick4], VElat1[tick5]), digits=2, nsmall=2), line=3, tick=FALSE)
    axis(1, at=RRtgrid, labels=format(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]), digits=2, nsmall=2), line=4.5, tick=FALSE)
  }
  mtext(expression(RR[t]), 1, line=1, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  mtext(expression(VE[0]^{lat}), 1, line=2.8, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  if(biomType == "Dichotomous") {
    mtext(expression(VE[2]^{lat}), 1, line=4.3, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  } else {
    mtext(expression(VE[1]^{lat}), 1, line=4.3, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
    mtext(expression(VE[2]^{lat}), 1, line=5.8, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  }
  mtext(expression("Vaccine Group CoR Relative Risk"~ RR[t]),1, line=7, at=mean(RRtgrid), cex=1.2)

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRt, power[,i], lty=i, col=colors[i], lwd=3)
  }

  abline(h=alpha/2,lty=3)

  legend(x="topright",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)

  title(main=bquote(paste("Power to Detect a "~.(biomType)~" CoR in Vaccine Recipients [2-sided ", alpha, "=", .(alpha), "]")))

}
