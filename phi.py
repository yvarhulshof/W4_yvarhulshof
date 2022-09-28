import numpy as np
from scipy.special import zetac, erf

def phi(mu,sigma):
	tau_s 		= 0.0005
	tau_r 		= 0.002
	tau_m 		= 0.01
	V_leak 		= -65
	V_reset 	= -65
	V_thresh 	= -50
	gamma = abs(zetac(0.5))/np.sqrt(2)
	s_tau = np.sqrt(tau_s/tau_m)
	N = len(mu)
	f = np.zeros(N)
	for n in range(1, N):
	    lower = (V_reset - V_leak - mu[n]) / sigma[n] + gamma * s_tau
	    upper = (V_thresh -V_leak - mu[n]) / sigma[n] + gamma * s_tau
	    x = np.linspace(lower,upper,100)
	    y = np.exp(x**2) * (1 + erf(x))
	    f[n] = (tau_r + tau_m * np.sqrt(np.pi) * np.trapz(x,y))**-1

	return f

