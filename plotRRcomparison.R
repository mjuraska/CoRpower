dataDir <- "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario5"
dataFile <- "ansScen5"
variableName <- "rho"
variableValues <- c(1,0.9,0.7,0.5)

plotRRcomparison <- function(dataDir, dataFile, variableName, variableValues, legendText=NULL) {
  load(paste0(file.path(dataDir[1], dataFile[1]),".RData"))
  power <- pwr$power
  RRt <- t(pwr$RRt)
  VElat0 <- pwr$VElat0
  VElat2 <- pwr$VElat2
  VElat1 <- rep(pwr$VEoverall, length(VElat0))
  alpha <- pwr$alpha
  RRlat2 <- pwr$RRlat2
  RRlat0 <- pwr$RRlat0
  ratio <- RRlat2/RRlat0
  
  ratio <- RRlat2/RRlat0
  
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
  
  if(is.null(legendText)) {
    legendText <- character(length(variableValues))
    for(i in 1:length(variableValues)) {
      legendText[i] <- paste0(variableName, " = ", variableValues[i])
    }
  }
  legend(x="topright",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)

}