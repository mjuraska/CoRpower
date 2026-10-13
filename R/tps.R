# Internal copy of the function 'tps' and its helpers from the R package 'osDesign' (version 1.8),
# Copyright (C) Sebastien Haneuse, Takumi Saegusa, Nilanjan Chatterjee, Norman Breslow,
# licensed under GPL (>= 3). Incorporated here to remove the dependency on 'osDesign', which is
# scheduled for archival on CRAN.
#
# Modified 2026-10-12 by Michal Juraska. Changes from the original:
#   - the estimation method is dispatched with switch() instead of eval(call())
#   - the internal helpers are renamed with the prefix 'tps' to avoid name clashes
#   - model.frame() is called as stats::model.frame()
#   - m$weight is spelled out as m$weights to avoid partial matching
#   - unused variables and commented-out code are removed, and the code is reformatted
# The numerical code is unchanged.
#
# References:
#   Breslow NE, Cain KC (1988). Logistic regression for two-stage case-control data. Biometrika 75: 11-20.
#   Breslow NE, Holubkov R (1997). Maximum likelihood estimation of logistic regression parameters
#     under two-phase, outcome-dependent sampling. JRSS-B 59: 447-461.
#   Haneuse S, Saegusa T, Lumley T (2011). osDesign: An R package for the analysis, evaluation, and
#     design of two-phase and case-control studies. Journal of Statistical Software 43(11): 1-29.

# 'tps' fits a logistic regression model to data arising from a two-phase (case-control or
# case-cohort) sampling design using the pseudo-likelihood (PL), weighted likelihood (WL), or
# maximum likelihood (ML) method
tps <- function(formula=formula(data), data=parent.frame(), nn0, nn1, group, contrasts=NULL,
                method="PL", cohort=TRUE, alpha=1){
  m <- match.call(expand.dots = FALSE)
  m$method <- m$contrasts <- m$nn0 <- m$nn1 <- m$group <- m$cohort <- m$alpha <- NULL
  m[[1]] <- quote(stats::model.frame)
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  X <- model.matrix(Terms, m, contrasts)

  if (length(nn0) != length(nn1)) stop("nn0 and nn1 should be of same length")
  if (length(nn0) != length(unique(group))) stop("Number of strata defined by group should be same as length of nn0")
  if (length(group) != nrow(X)) stop("Group and x are not compatible")
  if ((any(nn1 == 0)) || (any(nn0 == 0))) stop("Zero cell frequency at phase I")

  imeth <- charmatch(method, c("PL", "WL", "ML"), nomatch = 0)
  if (imeth == 0) stop("Method doesn't exist")

  if (is.matrix(Y) && (ncol(Y) > 1)){
    case <- as.vector(Y[, 1])
    N <- as.vector(Y[, 1] + Y[, 2])
  } else {
    case <- as.vector(Y)
    N <- rep(1, length(case))
  }

  fitFun <- switch(imeth, tpsPL, tpsWL, tpsML)
  out <- fitFun(nn0 = nn0, nn1 = nn1, x = X, N = N, case = case, group = group, cohort = cohort, alpha = alpha)
  names(out$coef) <- dimnames(X)[[2]]
  if (!is.null(out$cove)) dimnames(out$cove) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  if (!is.null(out$covm)) dimnames(out$covm) <- list(dimnames(X)[[2]], dimnames(X)[[2]])
  out$method <- method
  class(out) <- "tps"
  return(out)
}

# weighted likelihood (Horvitz-Thompson type) estimator
tpsWL <- function(nn0, nn1, x, N, case, group, cohort, alpha){
  nntot0 <- sum(nn0)
  nntot1 <- sum(nn1)
  u <- sort(unique(group))
  strt.id <- outer(group, u, FUN="==")
  strt.id <- matrix(as.numeric(strt.id), nrow=length(group), ncol=length(u))
  n1 <- apply(strt.id * case, 2, sum)
  n0 <- apply(strt.id * (N - case), 2, sum)
  if ((any(n0 == 0)) || (any(n1 == 0))) stop("Zero cell frequency at phase II")
  strt.id <- rbind(strt.id, strt.id)
  N <- c(case, N - case)
  case <- c(rep(1, length(case)), rep(0, length(case)))
  group <- c(group, group)
  x <- rbind(x, x)
  nstrta <- length(nn0)
  ofs <- 0
  if (!cohort) ofs <- log(sum(nn1)/sum(nn0)) - log(alpha)
  ofs <- rep(ofs, length(case))
  pw0 <- nn0[group]/n0[group]
  pw1 <- nn1[group]/n1[group]
  pw <- case * pw1 + (1 - case) * pw0
  Nw <- N * pw
  z <- rep(0, length(case))
  z[case == 1] <- Nw[case == 1]
  m <- glm(cbind(z, Nw - z) ~ -1 + x + offset(ofs), family = binomial, x=TRUE,
           control = glm.control(epsilon = 1e-06, maxit = 20))
  uj0 <- -t(m$x) %*% (strt.id * (1 - case) * m$fitted.values * N)
  uj1 <- t(m$x) %*% (strt.id * case * (1 - m$fitted.values) * N)
  g0 <- t(m$x * as.vector((1 - case) * m$fitted.values^2 * N * pw0^2)) %*% m$x
  g1 <- t(m$x * as.vector(case * (1 - m$fitted.values)^2 * N * pw1^2)) %*% m$x
  b0temp <- as.vector((nn0 * (nn0 - n0))/n0^3)
  b0 <- diag(b0temp, nrow=length(b0temp))
  b1temp <- as.vector((nn1 * (nn1 - n1))/n1^3)
  b1 <- diag(b1temp, nrow=length(b1temp))
  bb0temp <- as.vector((n0/(nntot0 * nn0)))
  bb0 <- diag(bb0temp, nrow=length(bb0temp))
  bb1temp <- as.vector((n1/(nntot1 * nn1)))
  bb1 <- diag(bb1temp, nrow=length(bb1temp))
  cov <- summary(m)$cov.unscaled
  cove <- g1 + g0 - uj0 %*% b0 %*% t(uj0) - uj1 %*% b1 %*% t(uj1)
  if (!cohort) cove <- cove - uj0 %*% bb0 %*% t(uj0) - uj1 %*% bb1 %*% t(uj1)
  cove <- cov %*% cove %*% cov
  return(list(coef = m$coef, cove = cove, fail = FALSE))
}

# pseudo-likelihood estimator of Breslow and Cain (1988)
tpsPL <- function(nn0, nn1, x, N, case, group, cohort, alpha){
  nntot0 <- sum(nn0)
  nntot1 <- sum(nn1)
  nstrata <- length(nn0)
  u <- sort(unique(group))
  strt.id <- outer(group, u, FUN="==")
  strt.id <- matrix(as.numeric(strt.id), nrow=length(group), ncol=length(u))
  n1 <- apply(strt.id * case, 2, sum)
  n0 <- apply(strt.id * (N - case), 2, sum)
  if ((any(n0 == 0)) || (any(n1 == 0))) stop("Zero cell frequency  at phase II")
  ofs <- log(n1/n0) - log(nn1/nn0)
  if (!cohort) ofs <- ofs + log(nntot1/nntot0) - log(alpha)
  ofs <- ofs[group]
  pw <- rep(1, length(case))
  m <- glm(cbind(case, N - case) ~ -1 + x + offset(ofs), family = binomial, weights = pw, x = TRUE,
           control = glm.control(epsilon = 1e-06, maxit = 20))
  fv <- m$fitted.values
  uj0 <- -t(m$x) %*% (strt.id * (N - case) * m$fitted.values)
  uj1 <- -t(m$x) %*% (strt.id * case * (1 - m$fitted.values))
  g0 <- t(m$x * (N - case) * fv^2) %*% m$x
  g1 <- t(m$x * case * (1 - fv)^2) %*% m$x
  a <- t(m$x) %*% (strt.id * (1 - fv) * fv * N)
  zz <- matrix(1, nrow = nstrata, ncol = nstrata)
  identity <- diag(rep(1, nstrata))
  b0 <- identity/as.vector(n0)
  b1 <- identity/as.vector(n1)
  bb0 <- identity/as.vector(nn0)
  bb1 <- identity/as.vector(nn1)
  if (!cohort){
    bb0 <- bb0 - (zz/nntot0)
    bb1 <- bb1 - (zz/nntot1)
  }
  aba <- a %*% bb0 %*% t(a) + a %*% bb1 %*% t(a)
  info <- t(m$x) %*% diag(m$weights) %*% m$x
  ghat <- info - a %*% b0 %*% t(a) - a %*% b1 %*% t(a)
  ghate <- g1 + g0 - uj0 %*% b0 %*% t(uj0) - uj1 %*% b1 %*% t(uj1)
  cov <- summary(m)$cov.unscaled
  # empirical variance-covariance matrix
  cove <- cov %*% (ghate + aba) %*% cov
  # model-based variance-covariance matrix
  covm <- cov %*% (ghat + aba) %*% cov
  return(list(coef = m$coef, covm = covm, cove = cove, fail = FALSE))
}

# maximum likelihood estimator of Breslow and Holubkov (1997) obtained by solving the concentrated
# Lagrangian equations with a modified Newton-Raphson algorithm
tpsML <- function(nn0, nn1, x, N, case, group, cohort, alpha, maxiter = 100){
  iter <- 0
  rerror <- 100
  nntot0 <- sum(nn0)
  nntot1 <- sum(nn1)
  nn <- nn0 + nn1
  nstrata <- length(nn0)
  nobs <- length(case)
  ncovs <- ncol(x)
  stpmax <- 1 * (ncovs + nstrata)
  u <- sort(unique(group))
  strt.id <- outer(group, u, FUN="==")
  strt.id <- matrix(as.numeric(strt.id), nrow=length(group), ncol=length(u))
  ee1 <- cbind(matrix(1, nstrata, 1), matrix(0, nstrata, nstrata))
  ee2 <- cbind(matrix(0, nobs, 1), strt.id)
  ee <- rbind(ee1, ee2)
  n1 <- apply(strt.id * case, 2, sum)
  n0 <- apply(strt.id * (N - case), 2, sum)
  if ((any(n0 == 0)) || (any(n1 == 0))) stop("Zero cell frequency at phase II")
  n1a <- c(nntot1, n1)
  n0a <- c(nntot0, n0)
  grpa <- c(rep(1, nstrata), (group + 1))
  na <- n0a + n1a
  yy <- c(nn1, case)
  identity <- diag(rep(1, nstrata))
  x0 <- cbind(identity, matrix(0, nrow = nstrata, ncol = ncovs))
  x <- cbind(-strt.id, x)
  xx <- rbind(x0, x)
  ofs <- log(n1a/n0a)
  ofs[1] <- ofs[1] - log(alpha)
  if (cohort) ofs[1] <- ofs[1] - log(nntot1/nntot0) + log(alpha)
  ofs <- ofs[grpa]
  repp <- c(nn1 + nn0, N)
  m <- glm(cbind(yy, repp - yy) ~ -1 + xx + offset(ofs), family = binomial, x = TRUE,
           control = glm.control(epsilon = 1e-10, maxit = 100))
  # initialize by Schill's estimates
  gamm0 <- gamm <- as.vector(m$coef)
  while ((iter < maxiter) && (rerror > 1e-10)){
    l <- tpsLagrange(gamm0, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy)
    h <- tpsLagrad(gamm0, nstrata, nobs, ncovs, nn1, nn0, n1a, n0a, ee, repp, grpa, n0, n1, xx, ofs)
    # Newton's direction
    p <- -solve(h) %*% l
    sp2 <- sqrt(sum(p^2))
    if (sp2 > stpmax) p <- p * (stpmax/sp2)
    gamm <- tpsLnsrch(gamm0, 0.5 * sum(l^2), t(h) %*% l, p, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy)
    rerror <- max(abs(gamm - gamm0)/pmax(abs(gamm0), 0.1))
    gamm0 <- gamm
    iter <- iter + 1
  }
  h <- tpsLagrad(gamm0, nstrata, nobs, ncovs, nn1, nn0, n1a, n0a, ee, repp, grpa, n0, n1, xx, ofs)
  fvv <- xx %*% gamm + ofs
  fvv <- as.vector(exp(fvv)/(1 + exp(fvv)))
  t1 <- nn1 - nn * fvv[1:nstrata]
  mu <- 1 - t1/n1
  t1 <- c(0, t1)
  mu <- c(1, mu)
  qn <- repp * n0a[grpa]
  qd1 <- n0a[grpa]
  qd2 <- (1 - mu[grpa]) * (n1a[grpa] - na[grpa] * fvv)
  q <- qn/(qd1 + qd2)
  gg <- t(xx) %*% ((1 - fvv) * fvv * q * ee)
  t00 <- (1 - fvv)^2 * q * ee
  d00 <- rep.int(1, nrow(t00)) %*% t00
  t01 <- (1 - fvv) * fvv * q * ee
  d01 <- rep.int(1, nrow(t01)) %*% t01
  t11 <- fvv * fvv * q * ee
  d11 <- rep.int(1, nrow(t11)) %*% t11
  d00 <- diag(as.vector(d00))
  d01 <- diag(as.vector(d01))
  d11 <- diag(as.vector(d11))
  tinfo <- (1 - fvv) * fvv * q * xx
  info <- t(xx) %*% tinfo
  cove <- solve(info)
  gig <- t(gg) %*% cove %*% gg
  zz1 <- cbind((gig + d00), (-gig + d01))
  zz2 <- cbind((-gig + d01), (gig + d11))
  zz <- rbind(zz1, zz2)
  zz <- solve(zz)
  gg <- cbind(gg, -gg)
  gig <- cove %*% gg %*% zz %*% t(gg) %*% cove
  # model-based variance-covariance matrix using the Aitchison and Silvey formula
  cov1 <- cove - gig
  # empirical variance-covariance matrix using within-group scores
  ta <- n0a * (n1a - t1)
  tb <- na * t1
  g <- (ta[grpa] * fvv)/(ta[grpa] + tb[grpa] * (1 - fvv))
  uj0 <- -t(m$x) %*% (ee * (repp - yy) * g)
  uj1 <- t(m$x) %*% (ee * yy * (1 - g))
  gg <- t(m$x * yy * (1 - g)^2) %*% m$x + t(m$x * (repp - yy) * g^2) %*% m$x
  gig <- gg - t(t(uj0)/n0a) %*% t(uj0) - t(t(uj1)/n1a) %*% t(uj1)
  h <- solve(h)
  cov2 <- t(h) %*% gig %*% h
  # adjustment for prospective first-stage sampling
  adjust <- (1/nntot0) + (1/nntot1)
  adjust <- matrix(adjust, nrow = (nstrata + 1), ncol = (nstrata + 1))
  if (cohort){
    cov1[1:(nstrata + 1), 1:(nstrata + 1)] <- cov1[1:(nstrata + 1), 1:(nstrata + 1)] + adjust
    cov2[1:(nstrata + 1), 1:(nstrata + 1)] <- cov2[1:(nstrata + 1), 1:(nstrata + 1)] + adjust
  }
  # extract the variance-covariance matrices for the regression coefficients
  cov1 <- cov1[(nstrata + 1):(nstrata + ncovs), (nstrata + 1):(nstrata + ncovs)]
  cov2 <- cov2[(nstrata + 1):(nstrata + ncovs), (nstrata + 1):(nstrata + ncovs)]
  delta <- gamm[1:nstrata]
  b <- gamm[(nstrata + 1):(nstrata + ncovs)]

  # validation of the constraints: h_ij from eq. (9) in Breslow and Holubkov (1997)
  xx2 <- xx[(nstrata + 1):nrow(xx), ]
  pred <- xx2 %*% gamm
  help <- rep(0, nstrata)
  for (i in 1:nstrata){
    help[i] <- -sum(xx2[, i])
  }
  n1_new <- rep(n1, help)
  n0_new <- rep(n0, help)
  p1 <- n1_new * exp(pred)/(n0_new + n1_new * exp(pred))
  pred2 <- xx[1:nstrata, ] %*% gamm
  if (cohort){
    P <- exp(pred2)/(1 + exp(pred2))
  } else {
    P <- nntot1 * exp(pred2)/(nntot0 + nntot1 * exp(pred2))
  }
  T <- nn1 - (nn1 + nn0) * P
  T_new <- rep(T, help)
  q <- N/(n0_new + n1_new) * n1_new * n0_new/(n0_new * n1_new + T_new * (n1_new - (n1_new + n0_new) * p1))
  p0 <- 1 - p1
  h1 <- rep(0, nstrata)
  h0 <- rep(0, nstrata)
  start <- 1
  for (i in 1:nstrata){
    end <- start + help[i] - 1
    sum_1 <- t(p1[start:end]) %*% q[start:end]
    sum_0 <- t(p0[start:end]) %*% q[start:end]
    h1[i] <- n1[i] - (n1[i] + n0[i]) * sum_1
    h0[i] <- n0[i] - (n1[i] + n0[i]) * sum_0
    start <- start + help[i]
  }
  fail <- any(round(h0, 10) != 0) || any(round(h1, 10) != 0)
  if (fail) cat("\n WARNING: ML estimates don't satisfy appropriate constraints\n\n")
  return(list(coef = b, covm = cov1, cove = cov2, delta = delta, mu = mu, fail = fail))
}

# Lagrangian equations
tpsLagrange <- function(gamm0, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy){
  fv <- xx %*% gamm0 + ofs
  fv <- exp(fv)
  fv <- fv/(1 + fv)
  r <- fv[1:nstrata]
  r <- nn1 - (nn1 + nn0) * r
  r <- c(0, r)
  r1 <- (n1a - r) * n0a
  r2 <- (n0a + n1a) * r
  g <- r1[grpa] * fv
  g <- g/(r1[grpa] + r2[grpa] * (1 - fv))
  g <- t(xx) %*% ((1 - g) * yy - g * (repp - yy))
  return(g)
}

# gradient of the Lagrangian equations
tpsLagrad <- function(gamm0, nstrata, nobs, ncovs, nn1, nn0, n1a, n0a, ee, repp, grpa, n0, n1, xx, ofs){
  fv <- xx %*% gamm0 + ofs
  fv <- exp(fv)
  fv <- fv/(1 + fv)
  xxx <- xx[(nstrata + 1):(nstrata + nobs), ]
  grp <- grpa[(nstrata + 1):(nstrata + nobs)] - 1
  rep <- repp[(nstrata + 1):(nstrata + nobs)]
  e <- ee[(nstrata + 1):(nstrata + nobs), 2:(nstrata + 1)]
  pp <- fv[1:nstrata]
  g <- fv[(nstrata + 1):(nstrata + nobs)]
  r <- nn1 - (nn1 + nn0) * pp
  r0 <- r * (n0 + n1)
  r1 <- n0 * n1 * (nn0 + nn1) * (n0 + n1) * pp * (1 - pp)
  r2 <- n0 * (n1 - r)
  d <- as.vector(r1[grp] * g * (1 - g) * rep/(r2[grp] + r0[grp] * (1 - g))^2)
  de <- d * e
  g <- -t(xxx) %*% de
  g <- cbind(g, matrix(0, (nstrata + ncovs), ncovs))
  r <- c(0, r)
  r1 <- n0a * (n1a - r)
  r2 <- (n0a + n1a) * r
  r0 <- n0a * n1a * (n1a - r) * (n0a + r)
  fv <- r0[grpa] * fv * (1 - fv) * repp/(r1[grpa] + r2[grpa] * (1 - fv))^2
  fv <- t(xx) %*% (tpsHdp(fv, xx))
  g <- g - fv
  return(g)
}

# horizontal direct product of two matrices
tpsHdp <- function(x, y){
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (nrow(x) != nrow(y)) stop("row dimensions are not same")
  z <- x[, rep(1:ncol(x), rep(ncol(y), ncol(x)))] * y[, rep(1:ncol(y), ncol(x))]
  return(z)
}

# line search for the modified Newton method used by 'tpsML'
tpsLnsrch <- function(gamm0, f0, g0, p, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy){
  scale <- 1.0
  rel.error <- max(abs(p)/pmax(gamm0, 0.1))
  if (rel.error < 1e-07){
    gamm <- gamm0
  } else {
    gamm <- gamm0 + p
    l <- tpsLagrange(gamm, nstrata, nn1, nn0, n1a, n0a, grpa, repp, xx, ofs, yy)
    f <- 0.5 * sum(l^2)
    slope <- sum(g0 * p)
    # backtrack if the function does not decrease sufficiently
    if (f > f0 + 0.0001 * slope){
      scale <- -slope/(2 * (f - f0 - slope))
      scale <- min(scale, 0.5)
      scale <- max(0.05, scale)
    }
    gamm <- gamm0 + scale * p
  }
  return(gamm)
}
