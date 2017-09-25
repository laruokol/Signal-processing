# f_noise
#
# A function for generating times series
# following a stochastic 1/f process.
#
# INPUT:   len   time series length
#          beta  desired slope for power series. Setting
#                beta = 0 gives unstructured 'white' noise,
#                while beta < 0 results in increasingly 
#                'red'-shifted series. A cautionary note:
#                when beta tends to -Inf, the series converges 
#                to one half of a sine wave. 
#          plot  should the time series be plotted or not?
#
# (c) Lasse Ruokolainen -- June 2015
##########################################################

f_noise = function(len,beta,plot=FALSE){
	u = (0:floor(len/2))
	S_f = (u^2)^(beta/2)
	S_f[is.infinite(S_f)] = 0
	phi = runif(length(S_f))
	z = Re(fft(S_f^0.5 * (cos(2*pi*phi)+(1i)*sin(2*pi*phi)),inv=T))
	if(plot==F){	
	}else{
		plot(z,type='l',
		main=paste('spectral slope (beta) =',as.character(beta)))
	}
	return(z)
}