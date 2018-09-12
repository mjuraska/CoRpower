dataDir <- c("T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario1/K1", "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario1/K2", "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario1/K3")
dataFile <- c("ansScen1", "ansScen1", "ansScen1")
variableName <- "controlCaseRatio"
variableValues <- c(5,3,1)

plotPowerTrinary <- function(dataDir, dataFile, variableName, variableValues, legendText=NULL) {
  load(paste0(file.path(dataDir[1], dataFile[1]),".RData"))
  power <- pwr$power
  RRt <- pwr$RRt
  VElat0 <- pwr$VElat0
  VElat2 <- pwr$VElat2
  VElat1 <- rep(pwr$VEoverall, length(VElat0))
  alpha <- pwr$alpha
  
  if(length(dataDir)>1) {
    for(i in 2:length(dataDir)) {
      load(paste0(file.path(dataDir[i], dataFile[i]),".RData"))
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
  tick1 <- which.min(abs(RRt-RRtgrid[1]))
  tick2 <- which.min(abs(RRt-RRtgrid[2]))
  tick3 <- which.min(abs(RRt-RRtgrid[3]))
  tick4 <- which.min(abs(RRt-RRtgrid[4]))
  tick5 <- which.min(abs(RRt-RRtgrid[5]))
  axis(1, at=RRtgrid, labels=round(c(VElat0[tick1], VElat0[tick2], VElat0[tick3], VElat0[tick4], VElat0[tick5]),2), line=1.5, tick=FALSE)
  axis(1, at=RRtgrid, labels=round(c(VElat1[tick1], VElat1[tick2], VElat1[tick3], VElat1[tick4], VElat1[tick5]),2), line=3, tick=FALSE)
  axis(1, at=RRtgrid, labels=round(c(VElat2[tick1], VElat2[tick2], VElat2[tick3], VElat2[tick4], VElat2[tick5]),2), line=4.5, tick=FALSE)
  mtext(expression(RR[t]), 1, line=1, at=-0.05, cex=1.2)
  mtext(expression(VE[0]^{lat}), 1, line=2.8, at=-0.05, cex=1.2)
  mtext(expression(VE[1]^{lat}), 1, line=4.3, at=-0.05, cex=1.2)
  mtext(expression(VE[2]^{lat}), 1, line=5.8, at=-0.05, cex=1.2)
  mtext(expression("Vaccine Group CoR Relative Risk"~ RR[t]),1, line=7, at=mean(RRtgrid), cex=1.2)
  
  colors <- c("blue","orange","forest green","black","red","purple","yellow", "pink")
  for(i in 1:ncol(power)){
    lines(RRt[,1], power[,i], lty=i, col=colors[i], lwd=3)
  }
  
  abline(h=alpha/2,lty=3)
  
  if(is.null(legendText)) {
    legendText <- character(length(variableValues))
    for(i in 1:length(variableValues)) {
      legendText[i] <- paste0(variableName, " = ", variableValues[i])
    }
  }
  
  legend(x="topright",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)
  
  title(bquote(paste("Power to Detect a Trichotomous CoR in Vaccine Recipients [2-sided ", alpha, "=", .(alpha), "]")))
  
}
