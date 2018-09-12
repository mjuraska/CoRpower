dataDir <- "T:/vaccine/StephanieWu/CoR Power Package In Progress/Scenario9"
dataFile <- "ansScen9"
variableName <- "rho"
variableValues <- c(1,0.9, 0.7, 0.5)


# dataDir and dataFile must be scalars
# rho = 1
plotVEcurveContinuous <- function(dataDir, dataFile, variableName, variableValues, legendText=NULL, risk0) {
  
  # Computes kernel of D(x, alphalat) (more information in Appendix B), which is the logit term in the zero-equation involving alphalat
  computeKernel <- function(x, alpha, nus, risk1latnu, sigma2obs){
    rho <- 1
    piece1 <- exp(alpha*(1 - x/nus[1]))*(risk1latnu^(x/nus[1]))
    piece2 <- (1-risk1latnu)^(x/nus[1]) + piece1
    piece3 <- dnorm(x/(sqrt(rho*sigma2obs)))
    kernel <- (piece1/piece2)*piece3
    return(kernel)
  }
  
  # Equation whose root gives alphalat. Labeled U(alphalat) in Appendix B and based on eqn 13 in the manuscript
  alphaLatEqn <- function(alpha, nus, risk1latnu, sigma2obs, VEoverall, PlatVElowest, risk0){
    logitterm <- integrate(computeKernel, lower=nus[1], upper=6, alpha=alpha, nus=nus, risk1latnu=risk1latnu, sigma2obs=sigma2obs)$value
    ans <- 1-VEoverall - (PlatVElowest*risk1latnu + logitterm)/risk0
    return(ans)
  }
  
  
  load(paste0(file.path(dataDir, dataFile),".RData"))
  power <- pwr$power
  RRc <- pwr$RRc
  alpha <- pwr$alpha
  sigma2obs <- pwr$sigma2obs
  PlatVElowest <- pwr$PlatVElowest
  VElowest <- pwr$VElowest
  VEoverall <- pwr$VEoverall
  
  rho <- 1
  nu <- sqrt(rho*sigma2obs)*qnorm(PlatVElowest)
  o <- length(VElowest)
  truebetas <- rep(NA,o)
  alphalatvect <- rep(NA,o)
  for (l in 1:o) {
    # find solutions alphalat and betalat by solving eqn (4) in Appendix B
    risk1latnu <- (1-VElowest[l])*risk0
    
    alphalatvect[l] <- uniroot(alphaLatEqn, lower=-10, upper=10, nus=nu, risk1latnu=risk1latnu, sigma2obs=sigma2obs, VEoverall=VEoverall, PlatVElowest=PlatVElowest, risk0=risk0)$root
    
    # solve for betalat:
    D <- risk1latnu
    truebetas[l] <- (log(D/(1-D)) - alphalatvect[l])/nu[1]
  }
  
  svect <- seq(-3,3,len=300)

  # Compute VE(s_1) vs s_1 for each fixed RRc value:
  inds <- round(seq(1,o,len=8))
  Fixedtruebetas <- truebetas[inds]
  FixedCoRRRs <- exp(Fixedtruebetas)
  
  for(i in 1:8) {
    linpart <- alphalatvect[inds[9-i]] + Fixedtruebetas[9-i]*svect
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
  
  tmp <- round(FixedCoRRRs[8],2)
  if (tmp!=1) { print("Warning for plot") }
  
  legend(x="bottomright",legend=c(as.expression(bquote("RR=1.00, "~VE[lowest]==.(round(VElowestvect[inds[8]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[7],2))~", "~VE[lowest]==.(round(VElowestvect[inds[7]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[6],2))~", "~VE[lowest]==.(round(VElowestvect[inds[6]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[5],2))~", "~VE[lowest]==.(round(VElowestvect[inds[5]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[4],2))~", "~VE[lowest]==.(round(VElowestvect[inds[4]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[3],2))~", "~VE[lowest]==.(round(VElowestvect[inds[3]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[2],2))~", "~VE[lowest]==.(round(VElowestvect[inds[2]],2)))),
                                  as.expression(bquote("RR="~.(round(FixedCoRRRs[1],2))~", "~VE[lowest]==.(round(VElowestvect[inds[1]],2))))),
         lty=c(1,2,3,4,5,6,7,8),col=c("brown","blue","dark grey","orange","red","purple","green","black"),lwd=2, cex=1.2)
  title("VE Curves for Different CoR Relative Risks per 1 SD in Vaccinees")
  text(2, 0.46, bquote(atop(VE[overall]==.(round(VEoverall,2))~", "~risk[0]==.(round(risk0,3)),
                            ~P[lowestVE]^{lat}==.(round(PlatVElowest,2)))),cex=1.3)
  
}

