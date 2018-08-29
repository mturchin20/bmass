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
	} else {
		c=matrix(0,nrow=1,ncol=p); 
		b=matrix(0,nrow=1,ncol=dd);
	}
	
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
	if(m==0){
		m = d-1
	}
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
