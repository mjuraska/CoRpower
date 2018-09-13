outDir <- "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario5"
outComputePower <- "ansScen5"
legendText <- paste0("rho = ", c(1, 0.9, 0.7, 0.5))  # must be specified in the same order as input parameters for computePower function
plotRRgradVE(outComputePower, outDir=outDir,legendText = legendText)

outComputePower <- list(pwr)
plotRRgradVE(outComputePower=outComputePower, legendText=legendText)

plotRRgradVE <- function(outComputePower, outDir=NULL, legendText) {
  if(is.list(outComputePower)) {
    pwr <- outComputePower[[1]] 
  } else if(is.character(outComputePower) & is.null(outDir)) {
    stop("outComputePower is a character vector so outDir needs to be specified")
  } else if(is.character(outComputePower)) {
    load(paste0(file.path(outDir[1], outComputePower[1]),".RData"))
  } else {
    stop("outComputePower must be of type list or character")
  }
  
  power <- pwr$power
  RRt <- t(pwr$RRt)
  alpha <- pwr$alpha
  RRlat2 <- pwr$RRlat2
  RRlat0 <- pwr$RRlat0
  ratio <- RRlat2/RRlat0
  
  if(length(outComputePower)>1) {
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
  
  legend(x="topleft",legend=legendText, lty=1:ncol(power),col=colors[1:ncol(power)],lwd=2,cex=1.2)

}
