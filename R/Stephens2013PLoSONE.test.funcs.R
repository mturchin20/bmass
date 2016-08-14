# Compute BFall (equations 12 and 13 in paper)
# for Y a d by N matrix, and X a N-vector of genotypes
# Added a parameter epsilon, with Psi = epsilon I, with the idea
# to check that if Y has exact linear dependencies then it does not effect the
# limiting BF
log10BFall = function(g,Y,sigmaa,N0=0,epsilon=0){
subset = complete.cases(t(Y)) & complete.cases(g)
Y=Y[,subset,drop=FALSE]
g=g[subset]
d = dim(Y)[1]
N = dim(Y)[2]
if(N0==0){N0 = d-1}
X = rbind(rep(1,N),g)
SXX = X %*% t(X)
  SXX[2,2] = SXX[2,2]+sigmaa^{-2}
  SYY = Y %*% t(Y)
  centerY = Y - apply(Y,1,mean)
  SYY0 = centerY %*% t(centerY)
  SYX = Y %*% t(X)
  LXX = chol(SXX)
  SYgX = SYY - SYX %*% chol2inv(LXX) %*% t(SYX)
  a = (d/2) * (-2*log(sigmaa) - 2*sum(log(diag(LXX)))+ log(N))
  b= (N+N0)/2 * (determinant(SYY0+diag(epsilon,d))$modulus - determinant(SYgX+diag(epsilon,d))$modulus)
return(list(lbf= (a+b)/log(10),a=a,b=b))
}



#Generate multivariate normal with given mean and covariance
#Generates a matrix with n columns, each column being d-variate N(mu,Sigma)
#[mu may be a matrix, in which case each column has a different mean]
rmvnorm = function(n,mu,Sigma){
  r = nrow(Sigma)
  Z = matrix(rnorm(n*r),nrow=r)
  L = chol(Sigma)
  return(mu + t(L) %*% Z)
}




#Generate random bivariate normal with unit variance and correlation rho
rbinorm = function(n,rho){
  return(rmvnorm(n,c(0,0),cbind(c(1,rho),c(rho,1))))
}

#note on prior
# Let di be the number of the d variables in each of the groups i=0,1,2
# we assume that the null model (all 0) has probability pi
# the sum of all alternatives (d1>=1) has probability 1-pi
# conditional on d1>=1, we assume that A=d1+d2 is uniform on 1...d,
# and conditional on A, d1 is uniform on 1 to A.
#conditional on d0, d1 and d2 we assume all d!/(d0! d1! d2!) configurations are equally likely
# thus Pr(given configuration) = pi if configuration is null,
# and (1-pi) (1/d) (1/(d1+d2)) d0! d1! d2!/d!
computeprior = function(z,pi){
  dvec = tabulate(z+1,nbin=3)
  d = length(z)
  return(ifelse(dvec[2]>0,(1-pi)*(1/d) * (1/(dvec[2]+dvec[3])) * 1/choose(d, dvec[1]) * 1/choose(d-dvec[1],dvec[2]), 0))
  
}

#picks out the partition gamma with all 1s (ie all in D)
allones = function(gamma){return(prod(gamma==1))}

#returns whether gamma matches a given value
gammamatch = function(gamma,g){return(identical(gamma,g))}


#This is one way to get VYY from Z scores
#Z is a p by d matrix;
#never used
#summaries.fromZscores = function(Z,n,f){
#  Zmax = apply(abs(Z),1,max)
#  nullsubset = (Zmax<2)
#  VYY = cor(Z[nullsubset,]) 
#}

#this function is similar to the "from summaries" function
#but allows us to deal with the fact that n may be different for each SNP
#note the "drop=FALSE" commands below stop R collapsing matrices into vectors inappropriately
#VYX \approx (1/n) Y'X is d by p
#VYY \approx (1/n) Y'Y is d by d
#VXX is a p-vector of the estimated variances of the SNP
logBF.fromVSummaries = function(VYX,VYY,VXX,U,D,n,m,d,sigmaa){

dd = sum(D)
du= sum(U)
p = dim(VYX)[2]

if(du>0){
LUU = chol(VYY[U,U,drop=FALSE]) # a du by du matrix
VUD = VYY[U,D,drop=FALSE]      #a du by dd matrix of the correlations of Yu with Yd

c = cbind(forwardsolve(t(LUU),VYX[U,,drop=FALSE]))#c solves LUU'c = phiU, c is a du by p matrix
b = cbind(forwardsolve(t(LUU), VUD))  # b is du by dd, and solves LUU' b = VUD, so b'b = VUD' LUU^-1 LUU'^-1 VUD = VUD' (LUU'LUU)^-1 VUD = VUD'VYYU^-1 VUD
} else{c=matrix(0,nrow=1,ncol=p); b=matrix(0,nrow=1,ncol=dd);}

C = VXX - colSums(c*c)

u = VYX[D,,drop=FALSE] - crossprod(b,c)

V0 = VYY[D,D,drop=FALSE] - crossprod(b)
L0 = chol(V0)
a = forwardsolve(t(L0),u)
lambda = sigmaa^(-2)  / (n*C)
k = as.numeric(1/(1+lambda))
return((dd/2) * log(1-k) - 0.5*(n+m-(d-sum(D)-sum(U)))*log(1-(k/C) *colSums(a*a)))
}

#computes the matrix of BFs for each value of sigmaa
#Z is a p by d matrix of Zscores
#n is a p vector of sample sizes
#f is a p vector of minor allele frequencies
#VYY is a d by d matrix of the estimated variance covariance matrix of Y
compute.allBFs.fromZscores = function(Z,VYY,n,f,sigmaa,pi0=0.5,m=0){
d = dim(Z)[2]
p = dim(Z)[1]
if(m==0){m = d-1}
VXX = 2*f *(1-f)
VYX = t(sqrt(VXX)*Z/sqrt(n))

prior = rep(0,3^d)
gamma=matrix(0,nrow=3^d,ncol=d)
lbf=list()
for(ss in 1:length(sigmaa)){
  lbf[[ss]] = matrix(0,nrow=3^d, ncol=p)
}

for(i in 0:(3^d-1)){
  for(j in 1:d){
    gamma[i+1,j]= (i %% 3^j) %/% 3^{j-1}
  }
  prior[i+1] = computeprior(gamma[i+1,],pi0)
  U = (gamma[i+1,]==0)
  D = (gamma[i+1,]==1)
  #print(U); print(D);
  if(prior[i+1]>0){
    for(ss in 1:length(sigmaa)){
      lbf[[ss]][i+1,] = logBF.fromVSummaries(VYX,VYY,VXX,U,D,n,m,d,sigmaa[ss])/log(10)
#note we just don't bother computing for models with prior = 0
    }
  } else {
    for(ss in 1:length(sigmaa)){
      lbf[[ss]][i+1,] = 0
    }
  }
}
prior[1] = pi0
return(list(lbf=lbf,gamma=gamma,prior=prior))
}

#Z is a p by d matrix of Zscores
#n is a p vector of sample sizes
#f is a p vector of minor allele frequencies
#VYY is a d by d matrix of the estimated variance covariance matrix of Y
#single-SNP BFs are computed for each SNP (row of Z)
#this particular funcion computes only the one BF Y_D | Y_U, G
logBF.rankone.Zmatrix.DgivenU = function(Z,VYY,n,f,Ucoord,Dcoord,sigmaa,pi0=0.5,m=0){
d = dim(Z)[2] #number of phenotypes
p = dim(Z)[1] #number of SNPs
if(m==0){m = d-1}
VXX = 2*f *(1-f)
VYX = t(sqrt(VXX)*Z/sqrt(n))

D=is.element(1:d,Dcoord) #true/false vector indicating which coordinates are in D
U = is.element(1:d,Ucoord)
return(logBF.fromVSummaries(VYX,VYY,VXX,U,D,n,m,d,sigmaa)/log(10))
}


#note the "drop=FALSE" commands below stop R collapsing matrices into vectors inappropriately
logBF.fromsummaries = function(phi,SYY,S11,U,D,n,m,d,sigmaa){
dd = sum(D)
du= sum(U)
p = dim(phi)[2]
if(du>0){
LUU = chol(SYY[U,U,drop=FALSE]) # a du by du matrix
phi0 = SYY[U,D,drop=FALSE]      #a du by dd matrix
c = cbind(forwardsolve(t(LUU),phi[U,,drop=FALSE]))#c solves LUU'c = phiU, c is a du by p matrix
b = cbind(forwardsolve(t(LUU), phi0))  # b is du by dd, and solves LUU' b = phi0, so b'b = phi0' LUU^-1 LUU'^-1 phi0' = phi0' (LUU'LUU)^-1 phi0 = phi0'SYYU^-1 phi0
} else{c=matrix(0,nrow=1,ncol=p); b=matrix(0,nrow=1,ncol=dd);}
C = S11 - colSums(c*c)
u = phi[D,,drop=FALSE] - crossprod(b,c)
RSS0 = SYY[D,D,drop=FALSE] - crossprod(b)
L0 = chol(RSS0)
a = forwardsolve(t(L0),u)
lambda = sigmaa^(-2)  / C
k = as.numeric(1/(1+lambda))
return((dd/2) * log(1-k) - 0.5*(n+m-(d-sum(D)-sum(U)))*log(1-(k/C) *colSums(a*a)))
}

#G is an n by p matrix of SNP genotypes
#Y is an n by d matrix of phenotypes
#single-SNP BFs are computed for each SNP (column of G)
#this particular funcion computes only the one BF Y_D | Y_U, G
logBF.rankone.matrix.DgivenU = function(G,Y,Ucoord,Dcoord,sigmaa,pi0=0.5,m=0){
if(is.null(dim(G))){G=cbind(G)} # turn a vector of genotypes into a matrix
subset = complete.cases(Y) & complete.cases(G)
Y=Y[subset,,drop=FALSE]
G=G[subset,,drop=FALSE]
n = dim(Y)[1]
d = dim(Y)[2]
p = dim(G)[2] #number of SNPs
if(m==0){m = d-1}
Y =scale(Y,center=T,scale=F)  #center Y and G to avoid needing intercept in regression
G = scale(G,center=T,scale=F)

phi = crossprod(Y,G) # this is t(Y) %*% G, a d by p matrix
SYY = crossprod(Y) # t(Y) %*% Y, a d by d matrix
S11 = colSums(G*G) # a p vector of ||g|| values
D=is.element(1:d,Dcoord) #true/false vector indicating which coordinates are in D
U = is.element(1:d,Ucoord)
return(logBF.fromsummaries(phi,SYY,S11,U,D,n,m,d,sigmaa)/log(10))
}


#G is an n by p matrix of SNP genotypes
#Y is an n by d matrix of phenotypes
#single-SNP BFs are computed for each SNP (column of G)
logBF.rankone.matrix = function(G,Y,sigmaa,pi0=0.5,m=0){
if(is.null(dim(G))){G=cbind(G)} # turn a vector of genotypes into a matrix
subset = complete.cases(Y) & complete.cases(G)
Y=Y[subset,,drop=FALSE]
G=G[subset,,drop=FALSE]
n = dim(Y)[1]
d = dim(Y)[2]
p = dim(G)[2] #number of SNPs
if(m==0){m = d-1}

Y =scale(Y,center=T,scale=F)  #center Y and G to avoid needing intercept in regression
G = scale(G,center=T,scale=F)

VYX = (1/n)*crossprod(Y,G) # this is (1/n) t(Y) %*% G, a d by p matrix
VYY = (1/n)*crossprod(Y) # (1/n) t(Y) %*% Y, a d by d matrix
VXX = (1/n)*colSums(G*G) # a p vector of (1/n) ||g|| values

prior = rep(0,3^d)
gamma=matrix(0,nrow=3^d,ncol=d)
lbf = matrix(0,nrow=3^d, ncol=p)

for(i in 0:(3^d-1)){
  for(j in 1:d){
    gamma[i+1,j]= (i %% 3^j) %/% 3^{j-1}
  }
  prior[i+1] = computeprior(gamma[i+1,],pi0)
  U = (gamma[i+1,]==0)
  D = (gamma[i+1,]==1)
  if(prior[i+1]>0){
    BF = 0
    for(ss in 1:length(sigmaa)){
      BF = BF+exp(logBF.fromVSummaries(VYX,VYY,VXX,U,D,n,m,d,sigmaa[ss]))#note we just don't bother computing for models with prior = 0
    }
    lbf[i+1,]=log(BF/length(sigmaa))
  } else {lbf[i+1,] = 0}
}
prior[1] = pi0

BF=exp(lbf)
posterior = prior[-1]*BF[-1,,drop=FALSE]
normalize=function(x){return(x/sum(x))}
posterior = apply(posterior,2,normalize)
p0=t(gamma[-1,]==0) %*% posterior
p1=t(gamma[-1,]==1) %*% posterior
p2=t(gamma[-1,]==2) %*% posterior

#divide by log(10) to convert everything to log base 10
lbfav = log(colSums(prior[-1]*exp(lbf[-1,,drop=FALSE]))/sum(prior[-1]))/log(10)
lbfuni.comps = apply(lbf[rowSums(gamma)==(2*d-1),,drop=FALSE],2,rev)/log(10)
lbfuni = log10(apply(10^lbfuni.comps,2,mean))
lbfall = lbf[which.max(apply(gamma,1,allones)),,drop=FALSE]/log(10)
lbf = lbf/log(10)
lbfuni.1 = lbfuni.comps[1] #these 4 are for the purposes of easily generating 2d simulation results in the paper
lbfuni.2 = lbfuni.comps[2]
gamma1given2 = rep(2,d)
gamma1given2[1]=1
gamma1given2[2]=0
lbf1given2 = lbf[apply(gamma,1,gammamatch,g=gamma1given2),,drop=FALSE]
gamma2given1 = rep(2,d)
gamma2given1[1]=0
gamma2given1[2]=1
lbf2given1 = lbf[apply(gamma,1,gammamatch,g=gamma2given1),,drop=FALSE]
return(list(prior=prior,gamma=gamma,lbf=lbf,lbfav = lbfav,lbfuni=lbfuni,lbfall=lbfall,lbfuni.comps=lbfuni.comps,lbfuni.1 = lbfuni.1 , lbfuni.2 = lbfuni.2,lbf1given2 =lbf1given2, lbf2given1 = lbf2given1,  p0=p0,p1=p1,p2=p2))
}


#test the first PC of Y for association with X
#Y is an n by d matrix
#X is an n-vector
prcomp.test=function(X,Y,sigmaa=0.4){
  w = eigen(cor(Y))$vectors[,1]
  z = Y %*% cbind(w)                    #this is Y projected on first PC
  return(logBF.rankone.matrix(X,z,sigmaa)) #test for association with X
}
