from numpy import *
from pylab import *
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
import sys
import os
# set tickers etc
rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Binding energies']})
title(r'{\bf Binding energies}', fontsize=20)     
# read in data from file
BEdata = loadtxt("liquid.dat")
l = len(BEdata[:,0])
x = BEdata[:,2]
y = BEdata[:,3]
z = BEdata[:,4]
axis([0,270,-1, 10.0])
xlabel(r'$A$',fontsize=20)
ylabel(r'$\mathrm{BE}$ [MeV]',fontsize=20)
plot(x, y ,'ro', x, z, 'bo', markersize=7)
# Save the figure in a separate file
savefig('beexpliquid.pdf', format='pdf')
# Draw the plot to screen
show()
    
