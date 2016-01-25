#!/usr/bin/env python

'''
Plot data files. Must have the same specific form as for example 2D_10pt_rs01.txt.
Plots only the result with the highest number of shells.
'''

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

#default
s_system = '2DHEG'

filename = sys.argv[1]
lines = np.loadtxt(filename, comments="#", unpack=False) #delimiter=""

sh = []
nw = []
ep = []
err = []

def shells2NumOrbs(shells):

	numsswf = np.zeros(len(shells))
	if s_system is '2DHEG':
		import shells2DHEG as shn
	elif s_system is '3DHEG':
		import shells3DHEG as shn

	numsswf = np.zeros(len(shells))
	for i in range(len(shells)):
		numsswf[i] = shn.states[shells[i]-1]

	return numsswf

#load columns in textfile
for i in range (1, len(lines[:,0]) ):
	if lines[i,0] > lines[i-1,0]: #plot only the last result for each nr of shells
		sh.append(int(lines[i-1,0]))
		#nw.append(int(lines[i-1,1]))
		ep.append(float(lines[i-1,2]))
		err.append(float(lines[i-1,3]))
sh.append(int(lines[i,0]))
#nw.append(int(lines[i,1]))
ep.append(float(lines[i,2]))
err.append(float(lines[i,3]))

plt.figure()
x = shells2NumOrbs(sh)
#x = x[::-1]
plt.errorbar(x,ep,yerr=err,fmt='k.')
#plt.yscale('log')
plt.show()#block=False)

