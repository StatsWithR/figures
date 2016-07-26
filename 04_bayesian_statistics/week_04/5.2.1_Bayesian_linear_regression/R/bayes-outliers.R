##
## Functions for Chaloner and Brant (1988) Biometrika 75:651-659 ##
## "A Bayesian Approach to Outlier Detection and Residual Analysis  ##
##

Bayes.outlier.prob <- function(lmobj, k=3) {
	e <- residuals(lmobj)
	h <- hatvalues(lmobj)
	Q <- qr.Q(lmobj$qr)
	alpha <- (lmobj$df.residual)/2
	rate <- (lmobj$df.residual*(summary(lmobj)$sigma)^2)/2 
	n <- length(e)
	pr <- rep(0,n)
	prjoint <- matrix(0,n,n)
	for (i in 1:n){
          pr[i] = integrate(outlier.prob,lower=0,upper=Inf,
                    ehat=e[i],hii=h[i],alpha=alpha,rate=rate,nsd=k)$value
        }
return(list(e=e,hat=h,prob.outlier=pr))
}

outlier.prob <- function(phi, ehat,hii,alpha,rate, nsd) {
	z1 <- (nsd - ehat*sqrt(phi))/sqrt(hii)
	z2 <- (- nsd - ehat*sqrt(phi))/sqrt(hii)
	pr.phi <- (1 - pnorm(z1) + pnorm(z2))*dgamma(phi,shape=alpha, rate=rate)
	return(pr.phi)}

bivoutlier.prob <- function(phi, ehati,ehatj,hii,hjj, rhoij, alpha,rate, nsd) {

        z1i = (nsd - ehati*sqrt(phi))/sqrt(hii)
	z2i = (- nsd - ehati*sqrt(phi))/sqrt(hii)
	z1j = (nsd - ehatj*sqrt(phi))/sqrt(hjj)
	z2j = (- nsd - ehatj*sqrt(phi))/sqrt(hjj)
        corr.neg = corr.pos =  diag(2)
        corr.pos[1,2] = corr.pos[2,1] = rhoij
        corr.neg[1,2] = corr.neg[2,1] = -rhoij
        B11 = apply(cbind(z1i,z1j), 1,
                    function(x){pmvnorm(lower=x, corr=corr.pos)} )
        B22 = apply(cbind(-z2i,-z2j), 1,
                    function(x){pmvnorm(lower=x, corr=corr.pos)})
        B12 = apply(cbind(z1i,-z2j), 1,
                    function(x){pmvnorm(lower=x, corr=corr.neg)})
        B21 = apply(cbind(-z2i,z1j), 1,
                    function(x){pmvnorm(lower=x, corr=corr.neg)})
        
        binorm = B11 + B22 + B12 + B21
        binorm[is.na(binorm)] = 0
#	binorm <-  (pmvnorm(-z1i,-z1j,corr=corr.pos) +
#                   pmvnorm(z2i,z2j,corr=corr.pos) - 
#	            pmvnorm(-z1i,z2j,corr=corr.neg) -
#                   pmvnorm(z2i,-z1j, corr=corr.neg))
	pr.phi <- binorm*dgamma(phi,shape=alpha, rate=rate)
	return(pr.phi)}


Bayes.outlier.prob.joint <- function(lmobj, k=3, joint=FALSE) {
	e <- residuals(lmobj)
	h <- hatvalues(lmobj)
	Q <- qr.Q(lmobj$qr)
	alpha <- (lmobj$df.residual)/2
	rate <- (lmobj$df.residual*(summary(lmobj)$sigma)^2)/2 
	n <- length(e)
	prjoint <- matrix(0,n,n)
	for (i in 1:n){
          j = 1
          while (j < i ) {
           corrij = sum(Q[i,]*Q[j,])/sqrt(h[i]*h[j]) 
           prjoint[i,j] = integrate(bivoutlier.prob,lower=0,upper=Inf,
                   ehati=e[i],ehatj=e[j],hii=h[i],hjj=h[j],
                   rhoij=corrij,alpha=alpha,rate=rate,nsd=k)$value
           prjoint[j,i] =  prjoint[i,j]
           j = j + 1
           print(c(i,j,prjoint[j,i]))
         }
        }
return(prjoint)
}

