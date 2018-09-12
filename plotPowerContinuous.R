dataDir <- "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario9"
dataFile <- "ansScen9"
variableName <- "rho"
variableValues <- c(1,0.9, 0.7, 0.5)


# if varying Plat2 and Plat0, must specify own legendText
plotPowerContinuous <- function(dataDir, dataFile, variableName, variableValues, legendText=NULL) {
  load(paste0(file.path(dataDir[1], dataFile[1]),".RData"))
  power <- pwr$power
  RRc <- pwr$RRc
  rho <- pwr$rho
  alpha <- pwr$alpha
  
  if(length(dataDir)>1) {
    for(i in 2:length(dataDir)) {
      load(paste0(file.path(dataDir[i], dataFile[i]),".RData"))
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
  
  if(is.null(legendText)) {
    legendText <- character(length(variableValues))
    for(i in 1:length(variableValues)) {
      legendText[i] <- paste0(variableName, " = ", variableValues[i])
    }
  }
  
  legend(x="topright", legend=legendText, lty=1:ncol(power), col=colors[1:ncol(power)], lwd=2, cex=1.2)
  
  title(bquote(paste("Power to Detect a Normally Distributed CoR in Vaccine Recipients [2-sided ", alpha, "=", .(alpha), "]")))
  
}
