##########################################
# PATH Rotavirus subunit efficacy trial

# For the rotavirus trial immune responses are measured at time tau = 3 months,
# but because events are only counted starting at month 3.5, we set tau=3.5.
# The end of follow-up for the endpoint is tmax = 24 months
# (notation from the paper).
tau <- 3.5
taumax <- 24

# 75% overall cumulative VE between Month 3.5 and 24
# (for the rotavirus trial VE is vs. a hypothetical placebo group):
VEoverall <- 0.75
RRoverall <- 1 - VEoverall

# With T the failure time, C the censoring time, X = min(T,C), and 
# Z vaccination status (vaccine or hypothetical placebo), we have
# RRoverall = P(tau < T <= taumax|Z=1)/P(tau < T <= taumax|Z=0) 

RRlatmed <- rep(RRoverall,100)
RRlatlo <- seq(1,RRoverall,len=100)

# Number enrolled in the subunit vaccine group:
Nv <- 4600

# Cumulative failure rate in the hypothetical placebo arm  
# between Month tau and taumax (which is rate0 in the paper):
#rate0 <- 0.034

# Cumulative failure rate in the hypothetical placebo arm
# between Month 0 and 24:
#rate0Mo0totaumax <- rate0*(1 + tau/taumax)
rate0Mo0totaumax <- 0.034

# Specify VE in the early period between enrollment and time tau --
# specify it half way between 0 and VEoverall:
VEoverall0totau <- VEoverall/2
RRoverall0totau <- 1 - VEoverall0totau  

# Cumulative dropout rate between Month 0 and 24:
Dropoutrate <- 0.10

# For Step 8, calculate the number of observed primary endpoint cases in a vaccine group 
# registered between Month 3.5 and 24, which is called ncases in the paper.  
# Calculate this under the following assumptions:
#
# 1. The failure time T and the censoring time C are independent.
# 2. T and C have exponential distributions in the hypothetical placebo arm
#    such that P(T <= t) = 1 - exp(-thetat t) and P(C <= t) = 1 - exp(-thetac) t).
# 3. RRoverall = P(tau < T <= t|Z=1)/P(tau < T <= t|Z=0) for all t between
#                tau and taumax (this will only approximately hold).
#

# Calculate thetat and thetac on a monthly time-scale:
thetat <- -log(1-rate0Mo0totaumax)/taumax
thetac <- -log(1-Dropoutrate)/taumax

# Calculate the number of vaccinees at-risk at Month 3.5:
# Logic: = Number enrolled * P(X > tau|Z=1)
#        =        Nv       * [1-RRoverall0totau*P(T <= tau|Z=0)]*P(C > tau)
#        =        Nv       * [1-RRoverall0totau*(1-exp(-thetat*tau))]*exp(-thetac*tau)
Nvtau <- Nv*(1-RRoverall0totau*pexp(tau,thetat))*(1-pexp(tau,thetac))

# where pexp(t,theta) is the cumulative distribution function of an exponential random
# variable with rate parameter theta.

# Number of observed rotavirus endpoints between Month 3.5 and 24 in vaccinees:  
# Logic: = Number vaccinees at-risk at Mo 3.5 * P(T<=24,T<=C|X>3.5,Z=1)
#        = Nvtau*\int_0^infty P(T<=min(c,taumax),C=c|X>tau,Z=1)dc
#        = Nvtau*\int_tau^infty P(3.5<T<=min(c,taumax),C=c|Z=1)dc / P(X>tau|Z=1)
#        = Nvtau*{\int_tau^infty P(3.5<T<=min(c,taumax)|Z=1)P(C=c)dc}/P(X>tau|Z=1).

# Now, the above integral in { } can be written as
#   \int_tau^infty P(3.5<T<=min(c,taumax)|Z=1)P(C=c)dc
# = \int_tau^taumax P(3.5<T<=c|Z=1)P(C=c)dc + P(3.5<T<=taumax|Z=1)\int_taumax^infty P(C=c)dc
# = RRoverall*\int_tau^taumax{[exp(-thetat*tau)-exp(-thetat*c)]
#            *thetac*exp(-thetac*c)dc} (term1)
# + RRoverall*[pexp(taumax,thetat)-pexp(tau,thetat)]*[1 - pexp(taumax,thetac)]. (term2)

# Term (1) is calculated to be:
#   RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
# - thetac*\int_tau^taumax exp(-(thetat+thetac)c)dc]
# = RRoverall*[exp(-thetat*tau)*[pexp(taumax,thetac)-pexp(tau,thetac)]
# - (thetac/(thetac+thetat))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))

# Lastly, the denominator (call it term3) equals:
# P(X>tau|Z=1) = [1-RRoverall0totau*P(T<=tau|Z=0)]P(C>tau)
#              = [1-RRoverall0totau*pexp(tau,thetat)][1-pexp(tau,thetac)].

# Putting it together
term1 <- RRoverall*(exp(-thetat*tau)*(pexp(taumax,thetac)-pexp(tau,thetac))
                    - (thetac/(thetac+thetat))*(pexp(taumax,thetac+thetat)-pexp(tau,thetac+thetat)))

term2 <- RRoverall*(pexp(taumax,thetat)-pexp(tau,thetat))*(1-pexp(taumax,thetac))
term3 <- (1-RRoverall0totau*pexp(tau,thetat))*(1-pexp(tau,thetac))

# ncases is notation in the paper:
ncases <- round(Nvtau*(term1+term2)/term3)
# Notation used in the R programs computerpower and computerpower2:
noobsatriskmotaucasesALL <- ncases

# Number of observed vaccinee controls eligible for potential sampling (those who
# reach the Month 24 visit without the rotavirus disease endpoint), which
# is calculated the same as for term3 except using taumax instead of tau
# (ncontrols is notation in the paper):
ncontrols <- round(Nv*(1-RRoverall*pexp(taumax,thetat))*(1-pexp(taumax,thetac)))

# Notation used in the R programs computerpower and computerpower2:
noobsatriskmotaucontrolsALL <- ncontrols

# Fraction of enrolled cases selected for the immunogenicity subset:
Fracimmunog <- 1.0  
noobsatriskmotaucasesPhase2 <- round(Fracimmunog*noobsatriskmotaucasesALL)

# Two-sided type I error rate of Wald tests:
alpha <- 0.05

# Control:case ratio:
controlcaseratio <- 5

# Parameters for the continuous marker case:
PlatVElowest <- 0.20 # For compatibility/plausibility, chosen to be smaller than the                                 # fraction not protected - 25%
VElowestvect <- seq(0,1-RRoverall,len=100)

# Measurement error parameters
sigma2obs <- 1
rhos <- c(1,0.9,0.7,0.5)

# M = number of Monte Carlo iterations for estimating power
M <- 1000

##############################################################
# Also consider 50% overall VE (for the live oral vaccine):
# VEoverall <- 0.50   
# etc.