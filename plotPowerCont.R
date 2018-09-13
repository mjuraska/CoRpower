#' @examples
#' 
#' outDir <- c("T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario9")
#' outComputePower <- c("ansScen9")
#' legendText <- paste0("rho = ", c(1, 0.9, 0.7, 0.5))  # must be specified in the same order as input parameters for computePower function
#' plotPowerCont(outComputePower, outDir=outDir,legendText = legendText)

#' outComputePower <- pwr
#' plotPowerCont(outComputePower=outComputePower, legendText=legendText)

plotPowerCont <- function(outComputePower, outDir=NULL, legendText) {
  
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
  RRc <- pwr$RRc
  rho <- pwr$rho
  alpha <- pwr$alpha
  
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
