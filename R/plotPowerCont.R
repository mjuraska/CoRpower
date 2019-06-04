#' Plotting of Power Curve versus Correlate of Risk Effect Size for Continuous Biomarkers
#'
#' Plots power (on the y-axis) to detect a correlate of risk effect size (on the x-axis) in the active treatment group for a continuous biomarker.
#' The correlate of risk effect size is quantified as the relative risk ratio of the clinical endpoint per standard deviation increase for a
#' noise-free biomarker.
#'
#' @param outComputePower either a list of lists containing output from \code{\link{computePower}} or a character vector specifying the \code{.RData} file(s) containing \code{\link{computePower}} output
#' @param outDir a character vector specifying path(s) to output \code{.RData} file(s), necessary if \cr
#' \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText a character vector specifying the entirety of the legend text. The order of the elements (i.e., parameter values) must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#' @param extendedLeg a logical value specifying if the extended legend with additional information about the control-to-case ratio, overall vaccine efficacy, number of cases, etc., is to be included. Default is \code{TRUE}.
#' @param xLegPos a number from \code{0} to \code{1} specifying the horizontal position of the extended legend, if applicable. A value of \code{0} produces text on the left side of the plot, \code{0.5} (default) produces text in the center, and \code{1} produces text on the right side.
#' @param yLegPos a number from \code{0} to \code{1} specifying the vertical position of the extended legend, if applicable. A value of \code{0} produces text at the bottom of the plot, \code{0.5} (default) produces text in the center, and \code{1} produces text at the top.
#' @param ySep a numeric value that specifies the spacing distance between lines in the extended legend, if applicable. Default is \code{0.7}.
#' @param margin a numeric vector of the form \code{c(bottom, left, top, right)}, which specifies the margins of the plot. Default is \code{c(7, 4, 3, 1)}. 
#'
#' @details The function's plot can be interpreted in conjunction with the output of \code{\link{plotVElatCont}} by
#' matching the CoR relative risk in the two plots and examining power compared to treatment (vaccine) efficacy.
#' This sheds light on the importance of overall vaccine efficacy on power and allows correlates of risk results
#' to be interpreted in terms of potential correlates of efficacy/protection.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' # Example scenario with continuous biomarker, where values of rho are varied
#'
#' # Set input parameters for computePower function
#' nCasesTx <- 10
#' nControlsTx <- 300
#' nCasesTxWithS <- 10
#' controlCaseRatio <- 5
#' VEoverall <- 0.75
#' risk0 <- 0.034
#' PlatVElowest <- 0.2
#' VElowest <- seq(0, VEoverall, len=5)
#' Plat0 <- P0 <- 0.2
#' Plat2 <- P2 <- 0.6
#' M <- 22
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.7, 0.4)
#' biomType <- "continuous"
#'
#' # Output from computePower function is stored in an object as a list of lists
#' pwr <- computePower(nCasesTx=nCasesTx, nCasesTxWithS=nCasesTxWithS, nControlsTx=nControlsTx,
#'                     controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                     PlatVElowest=PlatVElowest, VElowest=VElowest,
#'                     Plat0=Plat0, Plat2=Plat2, P0=P0, P2=P2, M=M, alpha=alpha,
#'                     sigma2obs=sigma2obs, rho=rho, biomType=biomType)
#'
#' # Set parameters for plotPowerCont function
#' # outComputePower is a list of lists containing output from the computePower function
#' outComputePower <- pwr
#' legendText <- paste0("rho = ", c(1, 0.7, 0.4))
#' plotPowerCont(outComputePower=outComputePower, legendText=legendText)
#'
#' \dontrun{
#' # Output from computePower function is saved in RData files
#' computePower(..., saveDir = "myDir", saveFile = "myFile.RData")

#' # outComputePower is a character string specifying the file containing the
#' # computePower output
#' # outDir is a character string specifying the outComputePower file directory
#' outComputePower <- paste0("myFile_rho_", c(1, 0.7, 0.4), ".RData")
#' outDir <- "~/myDir"
#' legendText <- paste0("rho = ", c(1, 0.7, 0.4))
#' plotPowerCont(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}}, \code{\link{plotVElatCont}}, \code{\link{plotPowerTri}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotPowerCont <- function(outComputePower, outDir=NULL, legendText,
                          extendedLeg=TRUE, xLegPos=0.5, yLegPos=0.5,
                          ySep=0.07, margin=c(7, 4, 3, 1)) {

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
  RRc <- pwr$RRc
  if (is.null(RRc)) {
    stop("Biomarker does not appear to be continuous. Consider using plotPowerTri() function for trichotomous biomarkers.")
  }
  rho <- pwr$rho
  alpha <- pwr$alpha

  if(length(outComputePower) > 1) {
    for(i in 2:length(outComputePower)) {
      if(is.list(outComputePower)) {
        pwr <- outComputePower[[i]]
      } else {
        if (length(outDir) != length(outComputePower)) {
          stop("outComputePower and outDir must have the same length")
        }
        load(paste0(file.path(outDir[i], outComputePower[i])))
      }
      addPower <- pwr$power
      power <- cbind(power, addPower)
    }
  }

  power <- as.matrix(power)

  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,mar=margin)
  plot(RRc,power[,1],ylim=c(0,1),type='n',xlab=expression("Vaccine Group CoR Relative Risk "~RR[c]~" per 1 SD Increase in S* with "~rho==1),ylab="Power",axes=FALSE)
  box()

  m <-length(RRc)
  RRcgrid <- seq(RRc[1],RRc[m],len=10)
  axis(1,at=RRcgrid,labels=round(RRcgrid,2))
  axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0),labels=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0))

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRc, power[,i], lty=i, col=colors[i], lwd=3)
  }

  abline(h=alpha/2,lty=3)

  legend(x="topright", legend=legendText, lty=1:ncol(power), col=colors[1:ncol(power)], lwd=2, cex=1.2)

  title(bquote(paste("Power to Detect a Normally Distributed CoR in Vaccine Recipients [2-sided ", alpha, "=", .(alpha), "]")))

  if (extendedLeg) {
    
    ### Add extra legend elements 
    label <- list()
    varyingArg <- pwr$varyingArg
    if (!(grepl("nCasesTx", varyingArg, fixed=TRUE))) {
      label$nCases <- bquote(n[cases]^S==~.(pwr$nCasesTxWithS))
    }
    if (!(grepl("controlCaseRatio", varyingArg, fixed=TRUE))) {
      label$controlCaseRatio <- bquote("controls:cases"==.(pwr$controlCaseRatio))
    }
    label$VEoverall <- bquote(VE[overall]==~.(round(pwr$VEoverall,2)))
    label$VElowest <- bquote(VE[lowest]==~.(round(max(pwr$VElowest), 2))~" to "~.(round(min(pwr$VElowest), 2)))
    if (!(grepl("PlatVElowest", varyingArg, fixed=TRUE))) {
      label$PlatVElowest <- bquote(P[lowestVE]^{lat}==~.(round(pwr$PlatVElowest, 2)))
    }
    if (!(grepl("rho", varyingArg, fixed=TRUE))) {
      label$rho <- bquote(rho==.(pwr$rho))
    }

    for (i in 1:length(label)) {
      text(min(RRcgrid) + xLegPos * (max(RRcgrid) - min(RRcgrid)), yLegPos - (i-1) * ySep, labels=label[[i]], pos=4)
    }
  }
  
}
