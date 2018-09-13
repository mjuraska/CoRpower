outDir <- "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario9"  
outComputePower <- "ansScen9"
risk0 <- 0.034
plotVElatContinuous(outComputePower, outDir=outDir, risk0 = risk0)

outComputePower <- list(pwr)
plotVElatContinuous(outComputePower = outComputePower, risk0 = risk0)

# outComputePower must be of length 1
# outDir must be of length 1
# rho = 1
# legend is determined by function
plotVElatContinuous <- function(outComputePower, outDir=NULL, risk0, legendText) {
  
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
  RRc <- pwr$RRc
  alpha <- pwr$alpha
  sigma2obs <- pwr$sigma2obs
  PlatVElowest <- pwr$PlatVElowest
  VElowest <- pwr$VElowest
  VEoverall <- pwr$VEoverall
  betaLat <- pwr$betaLat
  alphaLat <- pwr$alphaLat
  
  rho <- 1
  nu <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  o <- length(VElowest)
  
  svect <- seq(-3,3,len=300)

  # Compute VE(s_1) vs s_1 for each fixed RRc value:
  inds <- round(seq(1,o,len=8))
  fixedBetaLat <- betaLat[inds]
  fixedRRc <- RRc[inds]
  
  for(i in 1:8) {
    linpart <- alphaLat[inds[9-i]] + fixedBetaLat[9-i]*svect
    VEcurverr <- 1-(exp(linpart)/(1+exp(linpart)))/risk0
    VEcurverr[svect <= nu] <- VElowest[inds[9-i]]
    assign(paste0("VEcurverr",i), VEcurverr)
  }
  
  #ylims <- range(c(VEcurverr1,VEcurverr9,VEcurverr8,VEcurverr7,VEcurverr6,VEcurverr5,VEcurverr4,VEcurverr3))
  ylims <- c(-19,0.95)
  newylims <- log10(1-(ylims))
  
  
  par(cex.axis=1.4,cex.lab=1.4,cex.main=1.35,mar=c(4,5,4,4),oma=c(3,0,3,0),las=1)
  plot(svect,VEcurverr1,xlim=c(-3,3),ylim=c(-0.05,1),type='n',xlab=expression("True Biomarker X*=x* with "~rho==1~"in Vaccinees"), ylab=expression(VE[paste(x,"*")]^{lat}), axes=FALSE)
  axis(1)
  axis(2)
  
  colors <- c("brown","blue","dark grey","orange","red","purple","green","black")
  for(i in 1:8) {
    lines(svect, get(paste0("VEcurverr",i)), lty=i, col=colors[i], lwd=4)
  }
  
  tmp <- round(fixedRRc[8],2)
  if (tmp!=1) { print("Warning for plot") }
  
  legend(x="bottomright",legend=c(as.expression(bquote("RR=1.00, "~VE[lowest]==.(round(VElowest[inds[8]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[7],2))~", "~VE[lowest]==.(round(VElowest[inds[7]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[6],2))~", "~VE[lowest]==.(round(VElowest[inds[6]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[5],2))~", "~VE[lowest]==.(round(VElowest[inds[5]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[4],2))~", "~VE[lowest]==.(round(VElowest[inds[4]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[3],2))~", "~VE[lowest]==.(round(VElowest[inds[3]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[2],2))~", "~VE[lowest]==.(round(VElowest[inds[2]],2)))),
                                  as.expression(bquote("RR="~.(round(fixedRRc[1],2))~", "~VE[lowest]==.(round(VElowest[inds[1]],2))))),
         lty=c(1,2,3,4,5,6,7,8),col=c("brown","blue","dark grey","orange","red","purple","green","black"),lwd=2, cex=1.2)
  title("VE Curves for Different CoR Relative Risks per 1 SD in Vaccinees")
  text(2, 0.46, bquote(atop(VE[overall]==.(round(VEoverall,2))~", "~risk[0]==.(round(risk0,3)),
                            ~P[lowestVE]^{lat}==.(round(PlatVElowest,2)))),cex=1.3)
  
}

