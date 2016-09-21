
library(LearnBayes)
library(BayesFactor)
library(BAS)
source("BF-prop.R")

set.seed(8675309)
nsim=5000
sum = matrix(NA, nrow=nsim, ncol=6)
H = rbinom(nsim, p=.5, s=1)
na = nb = rep(100, nsim)
pa = rbeta(nsim, .5, .5)
pb = rbeta(nsim, .5, .5)
pc = rbeta(nsim, 1,1)
ya = rbinom(nsim, size=na, p=H*pa + (1- H)*pc)
yb = rbinom(nsim, size=nb, p=H*pb + (1- H)*pc)
for (i in 1:nsim) {

  # test.bayes.prop is BF of H0 to H1  so flip
  sum[i,1] = 1/bayes.prop.test(ya[i], na[i], yb[i], nb[i], alphaA=.5)$BF.H02H1
  #  LearnBayes gives BF of H1 to H0
  sum[i, 2] = exp(laplace(bfexch, 0, list(data=cbind(c(ya[i], yb[i]), c(na[i], nb[i])), K=1))$int)
  # intrinsic-moment is BF of H1 to H0
  sum[i, 3] = BF10.IM(ya[i], yb[i], na[i], nb[i])
  # Gunel-Dick is BF of H1 to H0
  sum[i, 4] = extractBF(contingencyTableBF(cbind(c(ya[i], yb[i]), c(na[i] - ya[i], nb[i] - yb[i])), 
                                 sampleType = "indepMulti", fixedMargin = "cols"))$bf

df.w = data.frame(y=c(1,1,0,0), x = c(1,0,1,0), w =c(ya[i], yb[i], na[i]- ya[i], nb[i]-yb[i]))
df = data.frame(y=c(rep(1, ya[i]),rep(1, yb[i]),rep(0, na[i]-ya[i]),rep(0, nb[i]-yb[i])),
                 x=c(rep(1, ya[i]),rep(0, yb[i]),rep(1, na[i]-ya[i]),rep(0, nb[i]-yb[i])))
                  
  #betaprior=CCH(a=2,b=2,s=0), 
  out.bas = bas.glm(y ~ x, data=df.w, weights=w, method="BAS",
                    betaprior=CCH(a=2, b=100, s=0),
                    modelprior=uniform(), family=binomial())
  sum[i,5] = exp(out.bas$logmarg[out.bas$size==2] - out.bas$logmarg[out.bas$size==1])
  sum[i,6] = anova(glm(y ~ x, data=df.w, weights=w, family=binomial()), test="Chisq")$Pr[2]

   }

ttt = data.frame(sum, H = factor(H))
confusion = table(ttt[,1] > 1, ttt$H)
print("default")
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion)

print("Laplace")
confusion = table(ttt[,2] > 1, ttt$H)
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion )

print("moment-intrinsic")
confusion = table(ttt[,3] > 1, ttt$H)
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion )

print("gunel-dickey")
confusion = table(ttt[,4] > 1, ttt$H)
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion) 

print("CCH")
confusion = table(ttt[,5] > 1, ttt$H)
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion )

print("pvalue")
confusion = table(ttt[,6] < .05, ttt$H)
(sum(confusion) - sum(diag(confusion)))/sum(confusion)
print(confusion )