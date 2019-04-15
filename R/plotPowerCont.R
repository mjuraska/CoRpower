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
#' M <- 22
#' alpha <- 0.05
#' sigma2obs <- 1
#' rho <- c(1, 0.7, 0.4)
#' biomType <- "continuous"
#'
#' # Output from computePower function is stored in an object as a list
#' pwr <- computePower(nCasesTx=nCasesTx, nCasesTxWithS=nCasesTxWithS, nControlsTx=nControlsTx,
#'                     controlCaseRatio=controlCaseRatio, risk0=risk0, VEoverall=VEoverall,
#'                     PlatVElowest=PlatVElowest, VElowest=VElowest, M=M, alpha=alpha,
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
#' outComputePower = c("myFile_rho_1.RData", "myFile_rho_0.7.RData", "myFile_rho_0.4.RData")
#' outDir = "~/myDir"
#' legendText <- paste0("rho = ", c(1, 0.7, 0.4))
#' plotPowerCont(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}}, \code{\link{plotVElatCont}}, \code{\link{plotPowerTri}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotPowerCont <- function(outComputePower, outDir=NULL, legendText) {

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
  rho <- pwr$rho
  alpha <- pwr$alpha

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

  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.2,las=1,oma=c(3,3,4,4))
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

}
