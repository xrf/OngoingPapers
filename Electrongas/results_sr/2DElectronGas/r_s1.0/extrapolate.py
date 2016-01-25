'''
1/M extrapolation

@numsswf 
@energies
@errors 
'''

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
from pylab import *
import sys
import os


class extrapolateEnergies:

	'''
	Constructor: input are numpy arrays of equal length
	@shells array with the number of energy shells in the basis
	@energies array with the number of energies
	@errors array with the errors of the energies
	@s_system physical system ('2DHEG' or '3DHEG')
	'''
	def __init__(self, shells, energies, errors, s_system):
		self.s_system = s_system
		self.numsswf = np.array([])
		self.energies = energies
		self.errors = errors
		self.popt = np.array([])
		self.pcov = np.array([])
		
		self.M2numsswf(shells)
		self.extrapolate()

	'''
	Reset arrays.
	'''
	def reset(self, shells, energies, errors, s_system):
		self.s_system = s_system
		self.numsswf = np.array([])
		self.energies = energies
		self.errors = errors
		self.popt = np.array([])
		self.pcov = np.array([])
		
		self.M2numsswf(shells)
		self.extrapolate()

	def M2numsswf(self, shells):
	
		self.numsswf = np.zeros(len(shells))
		if self.s_system is '2DHEG':
			import shells2DHEG as shn
		elif self.s_system is '3DHEG':
			import shells3DHEG as shn
		else:
			print('Error: value of s_system not valid.')
		
		for i in range(len(shells)):
			self.numsswf[i] = shn.states[shells[i]-1]

	'''
	Returns the energy assumet that it is on the form
	E(R) = a+b/x.
	'''
	def func(self,x,a,b):
		if self.s_system is '2DHEG':
			return a+b/np.float64(x)#**(3./2.))
		elif  self.s_system is '3DHEG':
			return a+b/np.float64(x)
	
	def extrapolate(self):
		it = 100000 #max iterations
		tol = 0.00001 #tolerance (see scipy manual)
		x = self.numsswf
		y = self.energies
		err = self.errors
		self.popt, self.pcov =\
				curve_fit(self.func, x, y, sigma=err, maxfev=it, gtol=tol)

	'''
	Find the extrapolated energy at infinity and the error 
	assumed that E(R,a,b) = a+b/M
	'''
	def eest(self):
		a = self.popt[0] 
		b = self.popt[1]

		p = np.matrix([[1, 0]]) #(d/da,d/db)(a+b/M) at M to infty
		err = np.sqrt(p*self.pcov*p.T)
		E = a
		print ('E: %f +/- %f')%(E,err)
		return E, err

	'''
	Plot the error relative to the extrapolated energy at infinity 
	and the extrapolated curve and energies with errorbars
	@lab plot legend
	@newfigure plot in new figure window
	'''
	def plottRelativeErrors(self, lab='', newfigure = False):	
		if newfigure:
			plt.figure()
		x = self.numsswf
		y = self.energies
		err = self.errors

		exp=-1.
		
		z = np.arange(x[0]/2,2*x[-1])
		res = self.popt[0]
		if lab is '':
			plt.plot(z**exp, - res+self.func(z, *self.popt))
		else:
			plt.plot(z**exp,self.func(z, *self.popt)-res, label='%s'%lab)
			plt.legend()
		plt.errorbar(x**exp,y-res,yerr=err,fmt='k.')
		
		plt.xscale('log')
		plt.yscale('log')

		tl = []
		for i in range (len(x)):
			tl.append('$%i^{-1}$'%x[i])
		plt.xticks(x**(exp), tl, rotation='-45')
		plt.xlabel('$1/M$')

		plt.xlim(1/(x[0]/2.), 1/(2.*x[-1]))
		plt.ylim( self.func(x[-1]*2., *self.popt)-res, self.func(x[0]/2., *self.popt)-res )
		
		plt.ylabel('$\Delta E$')

		plt.show(block=False)
		
#		plt.set_xlabel('$R$')
#		plt.set_xticklabels($1/%i^2/3$%x)#, rotation='45')
#		plt.set_xticks(x**exp)


'''
example

import numpy as np
y = np.array([-0.2346, -0.2486, -0.2523, -0.2533])
x = np.array([10,20,40,60])
err = np.array([.0003, .0002, .0005 , .0003])

# instantiate class

import extrapolate1OverM
obj = extrapolate1OverM.extrapolateEnergies(x,y,err,'2DHEG')

# error estimate
obj.eest()
#plott errors
obj.plottRelativeErrors('label', False)

# new vectors
x2=...
y2=...
err2=...
obj.reset(x2,y2,err2)

# An error estimate for the extrapolated energy are found on
# closed form (see masters thesis)
'''





dataSRG = loadtxt("plot.dat")
x = np.array(dataSRG[:,0], dtype=np.int8)
y=np.array(2*dataSRG[:,2])
err = np.ones(len(x))

obj = extrapolateEnergies(x,y,err,'2DHEG')

# error estimate
obj.eest()
#plott errors
obj.plottRelativeErrors('label', False)
