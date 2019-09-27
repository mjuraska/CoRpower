#' Plotting of Power Curve versus Correlate of Risk Effect Size for Continuous Biomarkers
#'
#' Plots power (on the y-axis) to detect a correlate of risk effect size (on the x-axis) in the active treatment group for a continuous biomarker.
#' The correlate of risk effect size is quantified as the odds ratio of the clinical endpoint comparing subgroups of active treatment recipients with a 1 standard deviation difference in a
#' noise-free biomarker response.
#'
#' @param outComputePower either a list of lists containing output from \code{\link{computePower}} or a character vector specifying the \code{.RData} file(s) containing \code{\link{computePower}} output
#' @param outDir a character vector specifying path(s) to output \code{.RData} file(s), necessary if \cr
#' \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText a character vector specifying the entirety of the legend text. The order of the elements (i.e., parameter values) must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#' @param legendTitle a character vector specifying the legend title if applicable (\code{NULL} by default)
#' @param extendedLeg a logical value specifying if the extended footnote legend with additional information about the control-to-case ratio, overall vaccine efficacy, number of cases, etc., is to be included. Default is \code{TRUE}.
#' @param verboseLeg a logical value specifying if the extended footnote legend shall use English words (\code{TRUE} by default) or mathematical notation used in Gilbert, Janes, and Huang (2016)
#' @param margin a numeric vector of the form \code{c(bottom, left, top, right)}, which specifies the margins of the plot. Default is \code{c(11, 7, 3, 1)}.
#'
#' @details If multiple levels are specified for the biomarker measurement error input argument \code{rho}, only the first level is used to determine
#' the \eqn{RR_c} values shown as x-axis tickmark labels.
#'
#' The function's plot can be interpreted in conjunction with the output of \code{\link{plotVElatCont}} by
#' matching the CoR relative risk in the two plots and examining power compared to treatment (vaccine) efficacy.
#' This sheds light on the importance of overall vaccine efficacy on power and allows correlates of risk results
#' to be interpreted in terms of potential correlates of efficacy/protection.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @references Gilbert P. B., Janes H., and Huang Y. (2016), Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials. \emph{Stat Med} 35(21):3745-59.
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
plotPowerCont <- function(outComputePower, outDir=NULL, legendText, legendTitle=NULL,
                          extendedLeg=TRUE, verboseLeg=TRUE, margin=c(11, 7, 3, 1)){

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
  plot(RRc,power[,1],ylim=c(0,1),type='n',xlab="",ylab="Power",axes=FALSE,cex.axis=1.3)
  xLabel <- ifelse(verboseLeg, "CoR Effect Size per 1 SD Increase in Noise-Free Biomarker", expression("CoR Effect Size"~RR[c]~"per 1 SD Increase in S* with"~rho==1))
  mtext(xLabel, side=1, line=2.5, cex=1.3)
  box()

  m <-length(RRc)
  RRcgrid <- seq(RRc[1],RRc[m],len=10)
  axis(1,at=RRcgrid,labels=round(RRcgrid,2))
  axis(2,at=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0),labels=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0))

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRc, power[,i], lty=i, col=colors[i], lwd=3.5)
  }

  abline(h=alpha/2,lty=3)

  legend(x="topright", legend=legendText, lty=1:ncol(power), col=colors[1:ncol(power)], lwd=3.5, cex=1.2, title=legendTitle)

  title(bquote(paste("Power to Detect a Continuous CoR in Vaccine Recipients [1-sided ", alpha, "=", .(alpha/2), "]")))

  if (extendedLeg) {

    ### Add extra legend elements
    label <- list()
    varyingArg <- pwr$varyingArg

    if (verboseLeg){

      idx <- 1
      label[[idx]] <- quote("Overall VE"==.(VEoverall)); idx <- idx + 1
      if (!(grepl("nCasesTx", varyingArg, fixed=TRUE))){
        if (!(grepl("controlCaseRatio", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote("Number controls"==.(pwr$controlCaseRatio * round(pwr$nCasesTxWithS))); idx <- idx + 1
        }
        label[[idx]] <- quote("Number cases"==.(round(pwr$nCasesTxWithS))); idx <- idx + 1
      }
      if (!(grepl("controlCaseRatio", varyingArg, fixed=TRUE))) {
        label[[idx]] <- quote("Controls:cases"==.(pwr$controlCaseRatio)*":1"); idx <- idx + 1
      }
      if (!(grepl("PlatVElowest", varyingArg, fixed=TRUE))){
        label[[idx]] <- quote("% lowest VE"==.(round(pwr$PlatVElowest * 100, 0))*"%"); idx <- idx + 1
      }
      if (!(grepl("rho", varyingArg, fixed=TRUE))) {
        label[[idx]] <- quote("Precision rho"==.(pwr$rho)); idx <- idx + 1
      }
      label[[idx]] <- quote("Observed var"==1); idx <- idx + 1

    } else {

      idx <- 1
      label[[idx]] <- quote("Overall VE"==.(VEoverall)); idx <- idx + 1
      if (!(grepl("nCasesTx", varyingArg, fixed=TRUE))){
        if (!(grepl("controlCaseRatio", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote(n[controls]^S==.(pwr$controlCaseRatio * round(pwr$nCasesTxWithS))); idx <- idx + 1
        }
        label[[idx]] <- quote(n[cases]^S==.(round(pwr$nCasesTxWithS))); idx <- idx + 1
      }
      if (!(grepl("controlCaseRatio", varyingArg, fixed=TRUE))) {
        label[[idx]] <- quote("Controls:cases"==.(pwr$controlCaseRatio)*":1"); idx <- idx + 1
      }
      if (!(grepl("PlatVElowest", varyingArg, fixed=TRUE))){
        label[[idx]] <- quote(P[lowestVE]^{lat}==.(round(pwr$PlatVElowest, 2))); idx <- idx + 1
      }
      if (!(grepl("rho", varyingArg, fixed=TRUE))) {
        label[[idx]] <- quote(rho==.(pwr$rho)); idx <- idx + 1
      }
      label[[idx]] <- quote(sigma[obs]^2==1); idx <- idx + 1
    }

    nLabels <- length(label)
    mtext(eval(substitute(bquote(lb1*";"~lb2*";"~lb3*";"~lb4), list(lb1=label[[1]], lb2=label[[2]], lb3=label[[3]], lb4=label[[4]]))), side=1, line=4.5, cex=1.1)
    if (nLabels==5){ mtext(eval(substitute(bquote(lb5), list(lb5=label[[5]]))), side=1, line=ifelse(verboseLeg, 5.5, 5.8), cex=1.1) }
    if (nLabels==6){ mtext(eval(substitute(bquote(lb5*";"~lb6), list(lb5=label[[5]], lb6=label[[6]]))), side=1, line=ifelse(verboseLeg, 5.5, 5.8), cex=1.1) }
    if (nLabels==7){ mtext(eval(substitute(bquote(lb5*";"~lb6*";"~lb7), list(lb5=label[[5]], lb6=label[[6]], lb7=label[[7]]))), side=1, line=ifelse(verboseLeg, 5.5, 5.8), cex=1.1) }
  }
}
