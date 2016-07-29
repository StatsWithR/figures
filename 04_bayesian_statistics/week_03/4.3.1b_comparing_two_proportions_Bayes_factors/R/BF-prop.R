# from Consonni, Forster, Rocca Statistical Science 2013
# default Bayes Factor
# note H1 is complex model prob not equal while H0 is that they are equal

BF10 = function(y1, y2, n1, n2, a=NULL,b=.5) {
  if (is.null(a)) {
    n = n1 + n2
    a11 = a12 = b*n1/n
    a21 = a22 = b*n2/n
  }
  else{ a11 = a[1,1]; a12=a[1,2]; a21=a[2,1]; a22=a[2,2]}
  logmarg1 = lbeta(y1 + a11, n1 - y1 + a12) + lbeta(y2 + a21, n2 - y2 + a22) -
    lbeta(a11, a12) - lbeta(a21, a22)
  logmarg0 =  lbeta(y1 + y2 + b, n1 + n2 - y1 - y2 + b) - lbeta(b, b)
  return(exp(logmarg1 - logmarg0))
}

#helper function 
Kah = function(a, h=0) {
  if (h == 0) {return(1)}
  else {
    a11 = a[1,1]; a12=a[1,2]; a21=a[2,1]; a22=a[2,2]
    j = 0:(2*h)
    Bratio = lbeta(a11 + j, a12) + lbeta(a21 + 2*h - j, a22) - lbeta(a11,a12) - lbeta(a21, a22)
    lfact = lgamma(2*h + 1) - lgamma(j+1) - lgamma(2*h - j + 1)
    neg1pow = (-1)^j
    return(sum(neg1pow*exp(Bratio + lfact)))
  }
}

# moment (non-local)  Bayes Factor
BF10.moment = function(y1,y2, n1, n2, a=NULL, b=.5, h=0) {
  if (is.null(a)) {
    a = matrix(rep(b,4), ncol=2)
    n = n1 + n2
    a[1,1] = a[1,2] = b*n1/n
    a[2,1] = a[2,2] = b*n2/n
  } 
  astar = a
  astar[1,1] = a[1,1] + y1 
  astar[1,2] = a[1,2] + n1 - y1
  astar[2,1] = a[2,1] + y2 
  astar[2,2] = a[2,2] + n2 - y2
  return(Kah(a=astar,h=h) * BF10(y1, y2, n1, n2, a=a, b=b)/Kah(a=a, h=h))
}

# marginal under H0
m0 = function(x1, x2, t1, t2, b) {
  return( exp(lbeta(b + x1 + x2, b + t1 + t2 -x1 -x2) - lbeta(b,b) + lgamma(t1+1) + lgamma(t2+1)-
                lgamma(x1+1) - lgamma(x2+1) - lgamma(t1 - x1 + 1) - lgamma(t2 - x2 + 1)))
}

# Balenced moment-based intrinsic prior
BF10.IM = function(y1, y2,n1, n2, a=NULL, b=.5, h=1, t1=4,t2=4) {
  if (is.null(a)) {
    a = matrix(rep(b,4), ncol=2)
    n = n1 + n2
    a[1,1] = a[1,2] = b*n1/n
    a[2,1] = a[2,2] = b*n2/n
  }
  BF = 0
  for  (x1 in 0:t1) {
    for (x2 in 0:t2) {
      x = matrix(c(x1, t1 - x1, x2, t2-x2), ncol=2, byrow=T)
      astar = a + x
      BF = BF + BF10.moment(y1, y2, n1, n2, a=astar, b=b, h=h)*m0(x1, x2, t1=t1, t2=t2, b=b)
    }
  }
  return(BF)
}

# default bayes factor with independent beta priors on theta_1, theta_2 and theta_0
# default settings as in Consonni et al

bayes.prop.test = function(yA, nA, yB, nB, pH0 = .5, alphaA=.25, betaA=NULL,alphaB=.25, betaB=NULL,
                           alphaC=NULL, betaC=NULL) {
  if (is.null(betaA)) betaA = alphaA
  if (is.null(betaB)) betaB = alphaB
  if (is.null(alphaC)) alphaC = alphaA + alphaB
  if (is.null(betaC))  betaC =  betaA + betaB
  logmarg0 = lbeta(yA + yB + alphaC, nA + nB - yA - yB + betaC) - lbeta(alphaC, betaC)
  logmarg1 = lbeta(yA + alphaA, nA - yA + betaA) + lbeta(yB + alphaB, nB - yB + betaB) -
    lbeta(alphaA, betaA) - lbeta(alphaB, betaB)
  BF0.1 = exp(logmarg0 - logmarg1)
  prior.oddsH1 = (1 - pH0)/pH0
  postprobH0 = 1/(1 + prior.oddsH1/BF0.1)
  return(list(postprobH0 = postprobH0, BF.H02H1 =  BF0.1))
}

# Example
#