#!/usr/bin/python
from numpy import *
from pylab import *
from matplotlib import rc, rcParams
import matplotlib.units as units
import matplotlib.ticker as ticker
import sys
import os


rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['New fits of nuclear forces']})

data6part = loadtxt("w_t_HO_COL_N6.dat")
data12part = loadtxt("w_t_HO_COL_N12.dat")
data20part = loadtxt("w_t_HO_COL_N20.dat")
data56part = loadtxt("w_t_HO_COL_N56.dat")

axis([0,1.0,0, 0.3])
xlabel(r'$\omega$',fontsize=20)
ylabel(r'$\langle T\rangle/\langle V\rangle$',fontsize=20)
plot(data6part[:,0],data6part[:,1]/(data6part[:,2]+data6part[:,3]) ,'b-*', data12part[:,0],  data12part[:,1]/(data12part[:,2]+data12part[:,3]),'r:.', data20part[:,0], data20part[:,1]/(data20part[:,2]+data20part[:,3]),'k--', data56part[:,0], data56part[:,1]/(data56part[:,2]+data56part[:,3]),'m:v', markersize=7)

#title(r'$N=6,\;\omega=0.5$', fontsize=20, horizontalalignment='center')
legend(('$N=6$','$N=12$','$N=20$','$N=56$'),
           'upper right', shadow=False, fancybox=False,prop={"size":18})
legend(loc='upper right')

#xticks( [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
# Save the figure in a separate file
savefig('qdots.pdf', format='pdf')
# Draw the plot to screen
show()
    
