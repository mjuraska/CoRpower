#' Plotting of Power versus Correlate of Risk Effect Size for Dichotomous and Trichotomous Biomarkers
#'
#' Plots power (on the y-axis) to detect a correlate of risk effect size (on the x-axis) in the active treatment group for a dichotomous or trichotomous biomarker. The correlate of risk effect size is quantified as
#' the relative risk of the clinical endpoint comparing subgroups of active treatment recipients with high and low biomarker response.
#'
#' @param outComputePower either a list of lists containing output from \code{\link{computePower}} or a character vector specifying the \code{.RData} file(s) containing \code{\link{computePower}} output
#' @param outDir a character vector specifying path(s) to output \code{.RData} file(s), necessary if \code{outComputePower} is a character vector. Default is \code{NULL}.
#' @param legendText a character vector specifying the entirety of the legend text. The order of the elements (i.e., parameter values) must match that of the \code{\link{computePower}} input parameters in order for legend labels to be accurate.
#' @param legendTitle a character vector specifying the legend title if applicable (\code{NULL} by default)
#' @param extendedLeg a logical value specifying if the extended footnote legend with additional information about the control-to-case ratio, overall vaccine efficacy, number of cases, etc., is to be included. Default is \code{TRUE}.
#' @param verboseLeg a logical value specifying if the extended footnote legend shall use English words (\code{TRUE} by default) or mathematical notation used in Gilbert, Janes, and Huang (2016)
#' @param margin a numeric vector of the form \code{c(bottom, left, top, right)}, which specifies the margins of the plot. Default is \code{c(11, 7, 3, 1)}.
#'
#' @details If multiple levels are specified for the biomarker measurement error input parameters (i.e., for \code{sens}/\code{spec} or \code{rho}) in \code{\link{computePower}}, only the first level is used to determine
#' the \eqn{RR_t} values shown as x-axis tickmark labels.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @references Gilbert P. B., Janes H., and Huang Y. (2016), Power/Sample Size Calculations for Assessing Correlates of Risk in Clinical Efficacy Trials. \emph{Stat Med} 35(21):3745-59.
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
#' outComputePower <- paste0("myFile_controlCaseRatio_", c(5, 3), ".RData")
#' outDir <- rep("~/myDir", 2)
#' legendText <- paste0("controls:cases = ", c("5:1","3:1"))
#' plotPowerTri(outComputePower, outDir=outDir, legendText = legendText)
#' }
#'
#' @seealso \code{\link{computePower}}, \code{\link{plotPowerCont}}
#'
#' @importFrom graphics abline axis box legend lines mtext par plot text title
#'
#' @export
plotPowerTri <- function(outComputePower, outDir=NULL, legendText, legendTitle=NULL,
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
  RRt <- pwr$RRt
  if (is.null(RRt)) {
    stop("Biomarker does not appear to be trichotomous. Consider using plotPowerCont() function for continuous biomarkers.")
  }
  VEoverall <- pwr$VEoverall
  VElat0 <- pwr$VElat0
  VElat2 <- pwr$VElat2
  alpha <- pwr$alpha
  Plat0 <- pwr$Plat0
  Plat2 <- pwr$Plat2
  Plat1 <- 1 - Plat0 - Plat2
  VElat1 <- (VEoverall - Plat0*VElat0 - Plat2*VElat2)/Plat1
  biomType <- ifelse(Plat1 == 0, "Binary", "3-Level")

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
  plot(RRt,power[,1],ylim=c(0,1),type='n',xlab="",ylab="Power",axes=FALSE, cex.axis=1.3)
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
  if(biomType == "Binary") {
    axis(1, at=RRtgrid, labels=format(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]), digits=2, nsmall=2), line=3, tick=FALSE)
  } else {
    axis(1, at=RRtgrid, labels=format(c(VElat1[tick1], VElat1[tick2], VElat1[tick3], VElat1[tick4], VElat1[tick5]), digits=2, nsmall=2), line=3, tick=FALSE)
    axis(1, at=RRtgrid, labels=format(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]), digits=2, nsmall=2), line=4.5, tick=FALSE)
  }
  mtext(ifelse(verboseLeg, "CoR RR", expression(RR[t])), 1, line=1, at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  mtext(ifelse(verboseLeg, "VE low", expression(VE[0]^{lat})), 1, line=ifelse(verboseLeg, 2.5, 2.8), at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  if(biomType == "Binary") {
    mtext(ifelse(verboseLeg, "VE high", expression(VE[2]^{lat})), 1, line=ifelse(verboseLeg, 4, 4.3), at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  } else {
    mtext(ifelse(verboseLeg, "VE med", expression(VE[1]^{lat})), 1, line=ifelse(verboseLeg, 4, 4.3), at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
    mtext(ifelse(verboseLeg, "VE high", expression(VE[2]^{lat})), 1, line=ifelse(verboseLeg, 5.5, 5.8), at=min(RRtgrid) - 0.15 * (max(RRtgrid) - min(RRtgrid)), cex=1.2)
  }

  varyingArg <- pwr$varyingArg

  # set the x-axis label
  # if rho is the varying input parameter, emphasize that RR_t is calculated for rho=1
  if (pwr$approach==2 & grepl("rho", varyingArg, fixed=TRUE)){
    if (verboseLeg){
      xLabel <- "CoR Effect Size [for Precision Rho = 1]"
    } else {
      xLabel <- expression("CoR Effect Size" ~ RR[t] ~~ "[for" ~ rho==1 * "]")
    }
  } else {
    if (verboseLeg){
      xLabel <- "CoR Effect Size"
    } else {
      xLabel <- expression("CoR Effect Size" ~ RR[t])
    }
  }
  mtext(xLabel, side=1, line=7, at=mean(RRtgrid), cex=1.3)

  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRt, power[,i], lty=i, col=colors[i], lwd=3.5)
  }

  abline(h=alpha/2,lty=3)

  legend(x="topright",legend=legendText, title=legendTitle, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=3.5,cex=1.2)

  title(main=bquote(paste("Power to Detect a "~.(biomType)~" CoR in Vaccine Recipients [1-sided ", alpha, "=", .(alpha/2), "]")))

  if (extendedLeg) {

    ### Add extra legend elements
    label <- list()

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
      if (!(grepl("Plat0", varyingArg, fixed=TRUE))){
        label[[idx]] <- quote("% low S = % low VE"==.(pwr$Plat0 * 100)*"%")
        label[[idx + 1]] <- quote("% high S = % high VE"==.(pwr$Plat2 * 100)*"%"); idx <- idx + 2
      }
      if (pwr$approach == 2){
        if (!(grepl("rho", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote("Precision rho"==.(pwr$rho)); idx <- idx + 1
        }
        label[[idx]] <- quote("Observed var"==1); idx <- idx + 1
      } else {
        if (!(grepl("sens", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote("Sens"==.(pwr$sens)~"; "~"Spec"==.(pwr$spec))
          label[[idx]] <- quote({FP^0==FN^2}==0)
        }
      }

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
      if (!(grepl("Plat0", varyingArg, fixed=TRUE))){
        label[[idx]] <- quote({P[0]^{lat}==P[0]}==.(pwr$Plat0))
        label[[idx + 1]] <- quote({P[2]^{lat}==P[2]}==.(pwr$Plat2)); idx <- idx + 2
      }
      if (pwr$approach == 2){
        if (!(grepl("rho", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote(rho==.(pwr$rho)); idx <- idx + 1
        }
        label[[idx]] <- quote(sigma[obs]^2==1); idx <- idx + 1
      } else {
        if (!(grepl("sens", varyingArg, fixed=TRUE))) {
          label[[idx]] <- quote("Sens"==.(pwr$sens)~"; "~"Spec"==.(pwr$spec))
          label[[idx]] <- quote({FP^0==FN^2}==0)
        }
      }
    }

    nLabels <- length(label)
    mtext(eval(substitute(bquote(lb1*";"~lb2*";"~lb3*";"~lb4), list(lb1=label[[1]], lb2=label[[2]], lb3=label[[3]], lb4=label[[4]]))), side=1, line=8.5, cex=1.1)
    if (nLabels==5){ mtext(eval(substitute(bquote(lb5), list(lb5=label[[5]]))), side=1, line=ifelse(verboseLeg, 9.5, 9.8), cex=1.1) }
    if (nLabels==6){ mtext(eval(substitute(bquote(lb5*";"~lb6), list(lb5=label[[5]], lb6=label[[6]]))), side=1, line=ifelse(verboseLeg, 9.5, 9.8), cex=1.1) }
    if (nLabels==7){ mtext(eval(substitute(bquote(lb5*";"~lb6*";"~lb7), list(lb5=label[[5]], lb6=label[[6]], lb7=label[[7]]))), side=1, line=ifelse(verboseLeg, 9.5, 9.8), cex=1.1) }
    if (nLabels==8){ mtext(eval(substitute(bquote(lb5*";"~lb6*";"~lb7*";"~lb8), list(lb5=label[[5]], lb6=label[[6]], lb7=label[[7]], lb8=label[[8]]))), side=1, line=ifelse(verboseLeg, 9.5, 9.8), cex=1.1) }

  }
}
