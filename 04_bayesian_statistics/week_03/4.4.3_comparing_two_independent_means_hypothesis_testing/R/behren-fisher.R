# Example Data
# ybar = c(52.10, 27.10)
# s2 = c(45.10, 26.40)^2
# n = c(22, 22)
# 
# BF = BFt.test(ybar, sd=sqrt(s2), n=n, sufficient=TRU, method="intrinsic", verbose=TRUE)
# 
# prior_H1 = 0.5
# prior_odds = prior_H1/(1 - prior_H1)
# post_odds = BF*prior_odds
# post_H1 = post_odds/(1 + post_odds)
# 
# 
# post1 = behren_fisher_intrinsic_BF_post_gibbs(ybar, s2, n, nsim=100000, burnin=10000)
# post2 = behren_fisher_intrinsic_BF_post(ybar, s2, n)


BFt.test = function(y, grp=NULL, n=NULL, sd.grp=NULL, prior.H1=.5, method="intrinsic", verbose=F,
                    nsim=10000, burnin=1000, thin=10,
                    plot=T, sufficient=F, jags=F, low.v=.0001, up.v=1000,m=4, c=2.3,
                    colHPD="gray", col.lines=NULL, ...) {

  if (jags) {
    require(R2jags)
  }
  else {
    require(cubature)
  }
  require(R2WinBUGS)
  
  if (is.null(grp) & !sufficient)  {
    stop("must supply grp if raw data rather than summary statistics are supplied")
  }
  if (sufficient & (is.null(n) | is.null(sd.grp))) {
    stop("if using sufficient statistics n and sd are required")
  }
  post.prob = 0

  if (!sufficient) {
    n <- length(y)
    if (!is.factor(grp))
      grp = as.factor(grp)
    grp.lab = levels(grp)
    means = tapply(y, grp, FUN = mean)
    sd = tapply(y, grp, FUN = sd)
    n = tapply(y, grp, FUN = length)
    
    data = list(
      ybarA = means[1],
      ybarB = means[2],
      nA = n[1],
      nB = n[2],
      ssA = (n[1] - 1) * sd[1] ^ 2,
      ssB = (n[2] - 1) * sd[2] ^ 2,
      eps = .000001
    )
  }
  else {
    means = y
    sd = sd.grp
    data = list(
      YbarA = y[1],
      YbarB = y[2],
      nA = n[1],
      nB = n[2],
      ssA = (n[1] - 1) * sd[1] ^ 2,
      ssB = (n[2] - 1) * sd[2] ^ 2,
      eps = .0001
    )
  }
if (verbose)  {
  print("means")
  print(means)
  print("standard deviations")
  print(sd)
  print("sample sizes")
  print(n)
}



  if (method == "independent.Jeffreys")  {
    muA = rt(nsim, df = (n[1] - 1)) * sd[1] / sqrt(n[1]) + means[1]
    muB = rt(nsim, df = (n[2] - 1)) * sd[2] / sqrt(n[2]) + means[2]
    diff = muA - muB
    diff.not0 = diff[diff != 0]
    ci = HPDinterval(as.mcmc(diff))
  }
  else if (method == "schwartz")  {
    BF = BFschwartz.H1.H0(means, sd ^ 2, n)
    post.prob = 1 - BF / (1 + BF)
    muA = rt(nsim, df = (n[1] - 1)) * sd[1] / sqrt(n[1]) + y[1]
    muB = rt(nsim, df = (n[2] - 1)) * sd[2] / sqrt(n[2]) + y[2]
    nzero = round(nsim * post.prob / (1 - post.prob))
    diff.not0 = muA - muB
    if (nzero > 0) {
      diff = c(muA - muB, rep(0, nzero))
    }
    else
      diff = diff.not0
    
    ci = HPDinterval(as.mcmc(diff))
  }
  
  else if (method == "matching") {
    model.M = function() {
      YbarA ~ dnorm(muA, nA * phiA)
      YbarB ~ dnorm(muB, nB * phiB)
      ssA ~ dgamma((nA - 1) / 2, phiA / 2)
      ssB ~ dgamma((nB - 1) / 2, phiB / 2)
      
      muA ~ dnorm(0, .0001)
      muB ~ dnorm(0, .0001)
      
      delta <- muA - muB
      
      phiA ~ dgamma(eps, eps)
      sigmaA <- pow(phiA, -.5)
      
      phiB <- pow(sigmaB, -2)
      # zero trick to define new prior distribution  logprior = -log(prior) + constant  > 0
      
      sigmaB ~ dunif(0, 1000)
      lambda <-
        10000 + 3 * log(sigmaB) - log(1 + pow(sigmaB / sigmaA, 2) * nA / nB)
      zero ~ dpois(lambda)
    }
    
    #M.model.file="BF-M.txt"
    #write.model(model.M, M.model.file)
    
    data$zero = 0
    params = c("muA", "muB", "delta", "sigmaA", "sigmaB")
    simM = jags(
      data = data,
      inits = NULL,
      parameters.to.save = params,
      model.file = model.M,
      n.iter = nsim,
      progress.bar = "none"
    )
    diff = simM$BUGSoutput$sims.matrix[, "delta"]
    diff.not0 = diff[diff != 0]
    ci = HPDinterval(as.mcmc(diff))
  }
  

  else if (method == "intrinsic" &  jags)  {
    model.Imeans = function()
    {
      muA <- mu * (1 - gamma) + muA1 * gamma
      muB <- mu * (1 - gamma) + muB1 * gamma
      precYA <- (phiA0 * (1 - gamma) + phiA1 * gamma)
      precYB <- (phiB0 * (1 - gamma) + phiB1 * gamma)
      
      
      YbarA ~ dnorm(muA, nA * precYA)
      YbarB ~ dnorm(muB, nB * precYB)
      ssA ~ dgamma((nA - 1) / 2, precYA / 2)
      ssB ~ dgamma((nB - 1) / 2, precYB / 2)
      
      muA1 ~ dnorm(mu, precA)
      muB1 ~ dnorm(mu, precB)
      
      precA <- 2 / (1 / phiA1 + 1 / phiA0)
      precB <- 2 / (1 / phiB1 + 1 / phiB0)
      
      phiA1 <- pow(sigmaA1, -2)
      phiB1 <- pow(sigmaB1, -2)
      
      sigmaA1 ~ dt(0, phiA0, 1) %_% T(0, )
      sigmaB1 ~ dt(0, phiB0, 1) %_% T(0, )
      
      phiA0 ~ dgamma(eps, eps)
      phiB0 ~ dgamma(eps, eps)
      
      v1 <- phiA1 / (phiA0)
      v2 <- phiB1 / (phiB0)
      
      mu ~ dnorm(0, .00001)
      gamma ~ dbern(prior.H1)
      
      delta <- muA - muB
    }
    
    
    #I.model.file="BF-I.txt"
    #write.model(model.Imeans, I.model.file)
    data$prior.H1 = prior.H1
    params = c("muA1",
               "muB1",
               "mu",
               "delta",
               "gamma",
               "sigmaA1",
               "sigmaB1",
               "v1",
               "v2")
    simM = jags(
      data = data,
      inits = NULL,
      parameters.to.save = params,
      model.file = model.Imeans,
      n.iter = nsim,
      progress.bar = "none"
    )
    diff = simM$BUGSoutput$sims.matrix[, "delta"]
    ci = HPDinterval(as.mcmc(diff))
    diff.not0 = diff[diff != 0]
    post.prob = mean(1 - simM$BUGSoutput$sims.matrix[, "gamma"])
    #  browser()
  }

  else if (method == "intrinsic" & !jags) {
    simM =  behren.fisher.post(
      ybar = means,
      s2 = sd ^ 2,
      n = n,
      low.v = low.v,
      up.v = up.v,
      m = m,
      prob.H1 = prior.H1,
      nsim = nsim,
      burnin = burnin,
      thin = thin,
      c = c
    )
    post.prob = 1 - simM$post.prob.H1
    nzero = round(nsim * post.prob / (1 - post.prob))
    diff.not0 = simM$mu1- simM$mu2
    if (nzero > 0) {
      diff = c(diff.not0, rep(0, nzero))
    }
    else {diff = diff.not0}
    
    ci = HPDinterval(as.mcmc(diff))
  }
    
else {warning("no method available; please check spelling"); return()}

#browser()
  if (plot) {
    if (is.null(col.lines)) {
      mydarkgrey = rgb(.5, .5, .5, name = "mydarkgrey", max = 1)
    }
    #   colHPD = rgb(86,155,189, name="myblue", max=256)
    
    
    dens = density(diff.not0)
    dens$y = (1 - post.prob) * dens$y / max(dens$y)
    #hist(diff.not0, prob=T, xlab=expression(mu[1] - mu[2]), main="posterior distribution", ...)
    plot(dens, ylim = c(0, 1), ...)
    polygon(c(ci[1], dens$x[dens$x > ci[1] & dens$x < ci[2]], ci[2]),
            c(0, dens$y[dens$x > ci[1] &
                          dens$x < ci[2]], 0),  col = colHPD)
    segments(0, 0, 0, post.prob, col = mydarkgrey)
  }
  
  return(list(
    mean.diff = mean(diff),
    ci = ci,
    diff = diff,
    post.prob.H0 = post.prob
  ))
}

logmarg.H0 = function(ybar, s2, n, low.mu=NULL, up.mu=NULL, m=4) {
# Marginal likelihood after integrating out phi1 and phi2 analytically
# This ignores some constants that will be constant across both hypotheses
# such as the gamma term in the sampling model for SS[i]
# ybar[i] ~ N(mu, (phi[i]n[i])^-1)
# ss[i] ~ Gamma((n[i]-1)/2, phi[i]/2
# p(mu, phi[1], phi[2]) proto 1/(phi[1]*phi[2]) 
#  returns the log marginal likelihood of H0

muhat = sum(ybar/s2)/sum(1/s2)
if (is.null(low.mu) | is.null(up.mu)) {
  offset = max(sqrt(s2))*m
  low.mu = muhat - offset
  up.mu = muhat + offset
}

    ss = s2*(n-1) 
    integrand = function(mu, ybar,ss,n, muhat) {
      
     ll = -0.5*(n[1]*log(0.5*(n[1]*((ybar[1] - muhat)^2) + ss[1])) + 
                n[2]*log(0.5*(n[2]*((ybar[2] - muhat)^2) + ss[2]))) +
          0.5*(log(n[1]) + log(n[2])) + lgamma(n[1]/2) + lgamma(n[2]/2)
     ll  = exp(-0.5*(n[1]*log(0.5*(n[1]*((ybar[1] - mu)^2) + ss[1])) + 
                     n[2]*log(0.5*(n[2]*((ybar[2] - mu)^2) + ss[2]))) +
                0.5*(log(n[1]) + log(n[2])) + lgamma(n[1]/2) + lgamma(n[2]/2) -
            ll)
          return(ll)
    } 
    log(integrate(integrand, low.mu, up.mu, ybar, ss, n, muhat)$value)
}


logmarg.H1 = function(ybar, s2, n, maxEval = 10 ^ 6,  low.v = .01,  up.v = 100, m = 4) {
  mu = (ybar[1] / s2[1] + ybar[2] / s2[2]) / (1 / s2[1] + 1 / s2[2])
  ss = s2 * (n - 1)
  offset = m * max(sqrt(s2))
  Lmu =  mu - offset
  Umu =  mu + offset
  log(
    adaptIntegrate(
      likelihood.mu.v1.v2,
      lowerLimit = c(Lmu, low.v, low.v),
      upperLimit = c(Umu, up.v, up.v),
      ybar,
      ss,
      n,
      maxEval = maxEval
    )$integral
  )
}



likelihood.mu.v1.v2 = function(x, ybar,ss,n ) {
  # Marginal likelihood after integrating out phi1 and phi2 analytically
  # This ignores some constants that will be constant across both hypotheses
  # such as the gamma term in the sampling model for SS[i]
  # ybar[i] ~ N(mu[i], (phi[i]n[i])^-1)
  # ss[i] ~ Gamma((n[i]-1)/2, phi[i]/2
  
  # mu[i] ~ N(mu, .5(sigma[i]^2 + tau[i]^2)
  # sigma[i] ~ C^+(0, tau[i]
  
  # p(mu, tau[1], tau[2]) proto 1/(tau[1]*tau[2])
  # returns the log marginal likelihood of H0
  
  
  mu = x[1]
  v1 = x[2]
  v2 = x[3]
  s2 = ss/(n - 1)
  muhat =(ybar[1]/s2[1] + ybar[2]/s2[2]) / (1/s2[1] + 1/s2[2])
  prec1 = 2 * n[1] / (n[1] + 2 + n[1] * v1)
  prec2 = 2 * n[2] / (n[2] + 2 + n[2] * v2)
  beta1 = .5 * (prec1 * (ybar[1] - mu) ^ 2 + ss[1])
  beta2 = .5 * (prec2 * (ybar[2] - mu) ^ 2 + ss[2])
  alpha1 = n[1]
  alpha2 = n[2]
  marg.loglike.H1 = -0.5*(n[1]*log(0.5*(n[1]*((ybar[1] - muhat)^2) + ss[1])) + 
                             n[2]*log(0.5*(n[2]*((ybar[2] - muhat)^2) + ss[2]))) +
    0.5*(log(n[1]) + log(n[2])) + lgamma(n[1]/2) + lgamma(n[2]/2)
  
  
  i = exp(
    .5 * (log(prec1) + log(prec2) - alpha1 * log(beta1) - alpha2 * log(beta2)) +
      lgamma(alpha1/2) + lgamma(alpha2/2) +
      -.5 * (log(v1) + log(v2)) - log(1 + v1) - log(1 + v2) - 2 * log(pi) - 
    marg.loglike.H1)
  
  return(i)
}


BFi.H1.H0 = function(ybar, s2, n, maxEval = 10^6, low.v = .001, up.v = 100, m = 4) {
  exp(logmarg.H1(ybar, s2, n, maxEval, low.v, up.v, m) - 
      logmarg.H0(ybar, s2, n, m = m))
}


behren.fisher.post.IS = function (ybar,
                                  s2,
                                  n,
                                  prob.H1 = .5,
                                  maxEval = 10 ^ 6,
                                  low.v = .001,
                                  up.v = 1000,
                                  m = 4,
                                  nsim = 10000) {
  odds = prob.H1 / (1 - prob.H1)
  BF = BFi.H1.H0(ybar, s2, n,  maxEval, low.v, up.v, m)
  post.odds = BF * odds
  post.prob.h1 = post.odds / (1 + post.odds)
  
  
  v1 = rf(nsim,  2, n[1])
  v2 = rf(nsim,  2, n[2])
  prec1 = 2 * n[1] / (2 + n[1] + n[1] * v1)
  prec2 = 2 * n[2] / (2 + n[2] + n[2] * v2)
  muhat = (ybar[1] / s2[1] + ybar[2] / s2[2]) / (1 / s2[1] + 1 / s2[2])
  sp = sqrt(s2[1] / n[1] + s2[2] / n[2])
  mu = rnorm(nsim, muhat, sd = sp)
  b1 = .5 * (s2[1] * (n[1] - 1) + (prec1 * n[1] * (ybar[1] - mu) ^ 2) / (prec1 + n[1]))
  b2 =  .5 * (s2[2] * (n[2] - 1) + (prec2 * n[2] * (ybar[2] - mu) ^ 2) / (prec2 + n[2]))
  phi1 = rgamma(nsim, (n[1] - 1) / 2, b1)
  phi2 = rgamma(nsim, (n[2] - 1) / 2, b2)
  mu1 = rnorm(nsim, (ybar[1] * n[1] + mu * prec1) / (n[1] + prec1),
              (phi1 * prec1) ^ {
                -1 / 2
              })
  mu2 = rnorm(nsim, (ybar[2] * n[2] + mu * prec2) / (n[2] + prec2),
              (phi2 * prec2) ^ {
                -1 / 2
              })
  
  logpost   = log(dnorm(ybar[1], mu1, 1 / sqrt((n[1] * phi1)))) +
    log(dnorm(ybar[2], mu2, 1 / sqrt((n[2] * phi2)))) +
    log(dnorm(mu1, mu,  sqrt(.5 * (1 + v1) / phi1))) +
    log(dnorm(mu2, mu,  sqrt(.5 * (1 + v2) / phi2))) +
    log(dgamma(s2[1] * (n[1] - 1), (n[1] - 1) / 2, phi1 / 2)) +
    log(dgamma(s2[2] * (n[2] - 1), (n[2] - 1) / 2, phi2 / 2)) +
    log(df(v1, 1 / 2, 1 / 2)) + log(df(v2, 1 / 2, 1 / 2)) +-log(phi1) - log(phi2)
  logprop  = log(dgamma(phi1, (n[1] - 1) / 2, b1)) +
    log(dgamma(phi2, (n[2] - 1) / 2, b2)) +
    log(df(v1, 2, n[1])) +
    log(df(v2, 2, n[2])) +
    log(dnorm(mu, muhat, sd = sp))
  log(dnorm(mu1, (ybar[1] * n[1] + mu * prec1) / (n[1] + prec1),
            1 / sqrt(phi1 * (n[1] + prec1))))
  log(dnorm(mu2, (ybar[2] * n[2] + mu * prec2) / (n[2] + prec2),
            1 / sqrt(phi1 * (n[2] + prec2))))
  wt = exp(logpost - logprop)
  wt = wt / sum(wt)
  
return(list(mu1=mu1, mu2=mu2, wt=wt, BF.H1.H0=BF, post.prob.H1=post.prob.h1))  
}

behren.fisher.post = function (ybar, s2, n, prob.H1 = .5, maxEval=10^6, 
                               low.v=.0001, up.v=1000, m=4,
                               nsim=10000, burnin=1000, thin=5, c=2.38) {
  odds = prob.H1 / (1 - prob.H1)
  #    browser()
  BF = BFi.H1.H0(ybar, s2, n,  maxEval, low.v, up.v, m)
  post.odds = BF * odds
  post.prob.h1 = post.odds / (1 + post.odds)
  
  # initialize
  theta.mat = matrix(NA, nrow = nsim, ncol = 7)
  colnames(theta.mat) = c("mu1", "mu2", "phi1", "phi2", "mu", "v1", "v2")
  
  theta = NULL
  theta$mu1 = ybar[1]
  theta$mu2 = ybar[2]
  theta$phi1 = 1 / s2[1]
  theta$phi2 = 1 / s2[2]
  theta$mu = .muhat(ybar, s2)
  theta$v1 = .01
  theta$v2 = .01
  
  for (i in 1:((nsim + burnin) * thin)) {
    # Gibbs  mu1, mu2, phi1, phi2 | v1, v2, mu
    theta = .update.mu.phi(theta, ybar, s2, n)
    #  theta = .update.mu(theta, ybar, s2, n, c)
    # Gibbs update for mu | mu1, mu2, phi1, phi2, v1, v2
    theta = .update.mu.gibbs(theta, ybar, s2, n)
    # random-walk Metropolis for log(v1)
    theta = .update.v1(theta, ybar, s2, n, c)
    # random-walk Metropolis  for log(v2)
    theta = .update.v2(theta, ybar, s2, n, c)
    # theta = .update.mu.v1.v2.marg(theta, ybar, s2, n, c)  NOT used now
    
    if (i > burnin & ((i - burnin) %% thin == 0)) {
      theta.mat[as.integer(i / thin) - burnin, ] = unlist(theta)
    }
  }
  return(list(mu1=theta.mat[,"mu1"], mu2=theta.mat[,"mu2"], BF.H1.H0=BF, 
              post.prob.H1=post.prob.h1, theta=theta.mat))  
}

.update.mu.v1.v2.marg = function(theta, ybar,s2,n,c) {
  
  sp = sqrt(sum(s2/n))
  sigma = c(sp*c, c,c)
  x = c(theta$mu, theta$v1, theta$v2)
  x.prop = c(x[1] + sp*c*rnorm(1), exp(log(x[2]) + c*rnorm(1)), 
            exp(log(x[3]) + c*rnorm(1)))

  R = log(likelihood.mu.v1.v2(x.prop, ybar, sqrt(s2*(n-1)), n)) -
      log(likelihood.mu.v1.v2(x, ybar, sqrt(s2*(n-1)), n)) + 
      sum(log(x.prop[2:3]) - log(x[2:3]))  
  if (R > log(runif(1))) {
    theta$mu = x.prop[1]
    theta$v1 = x.prop[2]
    theta$v2 = x.prop[3]
  }
  return(theta)
}
.update.mu.phi = function(theta, ybar,s2,n) {
  # block Gibbs update of mu1, mu2, phi1, phi2 | v1, v2, mu
  # draw phi1 | v1, mu  (marginal gamma)
  # draw phi2 | mid v2, mu  (marginal gamma)  
  # draw mu1 | phi1, v1, mu (conditional normal)
  # draw mu2 | phi2, v2, mu (conditional normal)
  
  v1 = theta$v1
  v2 = theta$v2
  mu = theta$mu
  prec1 =  2/(1 + v1)  #prior prec multiplier
  prec2 =  2/(1 + v2)
  b1 = .5*(s2[1]*(n[1] - 1) + (prec1*n[1]*(ybar[1] - mu)^2)/(prec1 + n[1]))
  b2 =  .5*(s2[2]*(n[2] - 1) + (prec2*n[2]*(ybar[2] - mu)^2)/(prec2 + n[2]))
  theta$phi1 = rgamma(1, (n[1]-1)/2, rate=b1)
  theta$phi2 = rgamma(1, (n[2]-1)/2, rate=b2)
  theta$mu1 = rnorm(1,(ybar[1]*n[1] + mu*prec1)/(n[1] + prec1),
                    sd=1/sqrt(theta$phi1*(n[1] + prec1)))
  theta$mu2 = rnorm(1,(ybar[2]*n[2] + mu*prec2)/(n[2] + prec2),
                    sd=1/sqrt(theta$phi2*(n[2] + prec2)))
  return(theta)
}
  
  
.muhat = function(ybar, s2)  {
 (ybar[1]/s2[1] + ybar[2]/s2[2])/(1/s2[1] + 1/s2[2])
}

.update.v1 = function(theta, ybar,s2,n, c) {
  theta.prop = theta
 theta.prop$v1 = exp(log(theta$v1) + c*rnorm(1))
  R = .logpost.Behren.Fisher.intrinsic(theta.prop,ybar,n,s2) -
      .logpost.Behren.Fisher.intrinsic(theta,ybar,n,s2) +
    log(theta.prop$v1) - log(theta$v1)
 
  if (R > log(runif(1))) {theta = theta.prop}
  return(theta)
}
.update.v2 = function(theta, ybar, s2, n, c) {
  theta.prop = theta
  theta.prop$v2 = exp( rnorm(1)*c + log(theta$v2))
  R = .logpost.Behren.Fisher.intrinsic(theta.prop,ybar,n,s2) -
      .logpost.Behren.Fisher.intrinsic(theta,ybar,n,s2) +
    log(theta.prop$v2) - log(theta$v2)
  if (R > log(runif(1))) { theta = theta.prop}
  return(theta)   
}
.update.mu = function(theta, ybar, s2, n, c) {
  theta.prop = theta
  sp = sqrt(sum(s2/n))
  theta.prop$mu = rnorm(1)*c*sp + theta$mu
  R = .logpost.Behren.Fisher.intrinsic(theta.prop,ybar,n,s2) -
      .logpost.Behren.Fisher.intrinsic(theta,ybar,n,s2) 
  if (R > log(runif(1))) { theta = theta.prop}
  return(theta)
}  

.update.mu.gibbs = function(theta, ybar, s2, n) {
prec = 2*(theta$phi1/(1 + theta$v1)  + theta$phi2/(1 + theta$v2))
mean = 2*(theta$mu1*theta$phi1/(1 + theta$v1) + theta$mu2*theta$phi2/(1 + theta$v2))/prec
theta$mu = rnorm(1, mean, sd=1/sqrt(prec)) 
  return(theta)
}

.logpost.Behren.Fisher.intrinsic = function(theta,ybar,n,s2)  {
  mu1 = theta$mu1
  mu2 = theta$mu2
  mu = theta$mu
  phi1= theta$phi1
  phi2 = theta$phi2
  v1 = theta$v1
  v2 = theta$v2
#  tau1 = sqrt(v1/phi1)
#  tau2 = sqrt(v2/phi2)
 logpost =  log(dnorm(ybar[1],mu1, 1/sqrt((n[1]*phi1)))) +
            log(dnorm(ybar[2],mu2, 1/sqrt((n[2]*phi2)))) +
            log(dnorm(mu1, mu,  sqrt(.5*(1 + v1)/phi1))) +
            log(dnorm(mu2, mu,  sqrt(.5*(1 + v2)/phi2))) +
            log(dgamma(s2[1]*(n[1]-1), (n[1] - 1)/2, phi1/2)) +
            log(dgamma(s2[2]*(n[2]-1), (n[2] - 1)/2, phi2/2)) +
           -.5*log(v1) - log(1 + v1) -.5*log(v2) - log(1 + v2) +
           - log(phi1) - log(phi2) 
 return(logpost)
}

# Schwartz approximation agrees for examples in paper
BFschwartz.H1.H0 = function (ybar, s2, n) {
mu =( ybar[1]/s2[1] + ybar[2]/s2[2])/(1/s2[1] + 1/s2[2])

exp(.5*(n[1]*log(1 + (ybar[1] - mu)^2/s2[1]) +
        n[2]*log(1 + (ybar[2] - mu)^2/s2[2]) -
        log(n[1]) - log(n[2]) + log(n[1] + n[2])))
}
