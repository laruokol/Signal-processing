### spatialNoise 
###
### A function to generate spatio-temporally autocorrelated
### stochastic variation, using a given covariance matrix. 
###
### INPUT:  C      spatial covariance matrix
###         Tmax   number of time points to be generated
###         sd     desired standard deviation (default: sd = 1)
###         ac     optional parameter setting the temporal
###                autocorrelation (default: ac = 0)
###         method should the series be generated using
###                reversed PCA ('pca', default) or Cholesky
###                factorisation ('chol')?
###
### (c) Lasse Ruokolainen -- December 2015
###
### Citation: Ruokolainen L (2013) PLoS ONE 8(8): e72325.
################################################################ 

spatialNoise = function(C,Tmax,sd = 1,ac = 0,method = 'pca'){

	#### EIGEN DECOMPOSITION OF C:
	# ----------------------------
	n = ncol(C)
	U = eigen(C)
		
	#### MAKE SURE MATRIX C IS OK:	
	# ----------------------------	
	# Check for symmetry of C:
	if(all(C == t(C)) == FALSE){
		stop('Covariance matrix is not symmetric')
	}
	# Check if C is positive definite:
	if(all(U$values > 0) == FALSE){
		stop('Covariance matrix is not positive definite')
	}

	#### GENERATE RANDOM SERIES:
	# --------------------------	
	X = matrix(rnorm(Tmax*n),Tmax,n)
    
   	#### GENERATE CORRELATED SERIES:
	# ------------------------------
	if(method == 'pca'){ # reversed PCA method
		Z = X %*% sqrt(diag(U$val)) %*% t(U$vec)
	}else{
		if(method == 'chol'){ # Cholesky factorisation method
			A = chol(C)
			Z = X %*% t(A)
		}else{
			stop('Non-supported method')
		}
	}
	
	#### TEMPORAL AUTOCORRELATION:
	# ----------------------------	
	if(ac != 0){
		W = matrix(0,Tmax,n)
		W[1,] = rnorm(n)
		for(t in 2:Tmax){
			W[t,] = ac*W[t-1,] + Z[t,]
		}
		Z = W
	}
	
	#### STANDARDISE VARIANCES:
	# -------------------------	
	Z = sd * apply(Z,2,function(x){(x-mean(x))/sd(x)})
	
	#### RETURN CORRELATED SERIES: 
	# ----------------------------	
	return(Z)
}