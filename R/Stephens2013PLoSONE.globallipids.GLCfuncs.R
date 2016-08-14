
#function to return the top hits exceeding some threshold,
#removing duplicates within 1Mb in a greedy way
indephits=function(score,chr,pos,T=5e5){
R=rep(0,length(score))
for(i in 1:length(score)){
  R[i] = score[i]==max(score[chr==chr[i] & pos<(pos[i]+T) & pos>(pos[i]-T)])
}
return(R)
}

normalize=function(x){return(x/sum(x))}

#remove max from each column to avoid overflow
centered.lbf = function(lbf){
  return(t(t(lbf) - apply(lbf,2,max)))
}


#compute posterior from prior and matrix of log10 bfs
#lbf is nmodels by p snps
#prior is nmodel vector
posteriorprob = function(lbf, prior){
lbf =centered.lbf(lbf) #remove max from each column to avoid overflow
post = apply(prior * 10^lbf,2,normalize)
return(post)
}

#note it is important that priorinit is 0 for models that are impossible
#since lbf is typically meaningless for those in our current implementation
em.priorprobs = function(lbf,priorinit,niter=10){
  pp = priorinit
  for(i in 1:niter){
    pp = apply(posteriorprob(lbf,pp),1,mean)
  }
  return(pp)
}

#compute marginal probabilities for being in U, D, I classes for each phenotype
#gamma is nmodel by ntraits (d)
#postprobs is an nsigma*nmodel by nsnps matrix composed from stacking nsigma matrices (one for each value of sigma) each nmodel by nsnps 

marginal.postprobs=function (postprobs, gamma, nsigma) 
{
  nmodel = dim(gamma)[1]
  d = dim(gamma)[2]
  gammarep = matrix(0, nrow = nsigma * nmodel, ncol = d)
  for (i in 1:d) {
    gammarep[, i] = rep(gamma[, i], nsigma)
  }
  return(list(pU = t(gammarep == 0) %*% postprobs, pD = t(gammarep == 
                                                            1) %*% postprobs, pI = t(gammarep == 2) %*% postprobs))
}



#compute average BF, averaged over prior
#note prior should just be on the non-null models here
lbf.av=function(lbf,prior){
  max = apply(lbf,2,max) #remove the max, and store it, to prevent overflow
  lbf=centered.lbf(lbf)
  bfav = prior * 10^lbf
  return(log10(apply(bfav,2,sum))+max)
}

# compute simple average of univariate bfs
# lbf is nmodel*nsigma by nsnps
#if multiple sigma values used, we just stack them on top of one another
#in lbf
lbf.uni=function(lbf,gamma){
  nsigma= dim(lbf)[1]/dim(gamma)[1]
  d = dim(gamma)[2]
  s = apply(gamma,1,sum)
  uni = (s == 2*d-1)
  max = apply(lbf,2,max)
  lbf = centered.lbf(lbf)
  bfuni = apply(10^lbf*uni,2, sum)/(d*nsigma) 
  return(log10(bfuni)+max)
}

# compute simple average of univariate bfs
# lbf is nmodel by nsnps
lbf.all=function(lbf,gamma){
  nsigma= dim(lbf)[1]/dim(gamma)[1]
  all = apply(gamma,1,allones)
  max = apply(lbf,2,max)
  lbf = centered.lbf(lbf)
  bfall = apply(10^lbf*all,2, sum)/nsigma 
  return(log10(bfall)+max)
}

