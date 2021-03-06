from pylab import *
import sys
import os
import numpy as np
from scipy.optimize import curve_fit

def func(x, a, b, c):
     return a+ b*(x**(-c))

filename = sys.argv[1]
startR 	 = int(sys.argv[2])

dataSRG = loadtxt(filename)


states=[ 2, 10, 18, 26, 42, 50, 58, 74, 90, 98, 114, 122, 138, 162, 178, 194, 202, 218, 226, 242, 258, 274, 290, 298, 322, 338, 354, 370, 386, 394, 426, 442, 450, 466, 482, 498, 506, 522, 554, 570, 586, 602, 610, 634, 650, 666, 682, 698, 714, 730, 746, 754, 770, 802, 810, 842, 858, 874, 882, 914, 930, 946, 962, 978, 994, 1010, 1018, 1034, 1058, 1090, 1106, 1122, 1138, 1154, 1186, 1202, 1218, 1226, 1242, 1266, 1282, 1314, 1330, 1346, 1362, 1394, 1418, 1434, 1450, 1466, 1482, 1498, 1514, 1522, 1538, 1554, 1586, 1594, 1610, 1642, 1658, 1690, 1706, 1722, 1738, 1754, 1770, 1778, 1802, 1834, 1850, 1866, 1882, 1898, 1930, 1946, 1962, 1978, 1994, 2010, 2018, 2066, 2082, 2098, 2114, 2138, 2170, 2186, 2202, 2218, 2234, 2250, 2258, 2274, 2306, 2322, 2354, 2370, 2402, 2418, 2434, 2450, 2458, 2474, 2490, 2514, 2530, 2546, 2562, 2578, 2610, 2626, 2642, 2658, 2706, 2722, 2738, 2746, 2778, 2810, 2826, 2850, 2866, 2882, 2898, 2914, 2930, 2946, 2962, 2978, 3010, 3026, 3034, 3066, 3082, 3098, 3130, 3162, 3194, 3210, 3218, 3234, 3266, 3282, 3298, 3306, 3338, 3370, 3386, 3402, 3418, 3450, 3466, 3482, 3498, 3514, 3530, 3562, 3578, 3586, 3602, 3626, 3658, 3674, 3706, 3722, 3738, 3754, 3770, 3786, 3802, 3834, 3850, 3866, 3882, 3922, 3938, 3954, 3986, 4002, 4018, 4034]

l = len(dataSRG)

for i in range(len(dataSRG)):
	if(dataSRG[i,0] == startR):
		l -= i
		startIndex = i
		break

x = ones(l)
y = ones(l)

for i in range(l):
	x[i] = states[(int)(dataSRG[i+startIndex,0]-1+0.1)]
	y[i] = dataSRG[i+startIndex,2]

print x

popt, pcov = curve_fit(func, x, y, sigma=None, maxfev = 200000, gtol=0.00001)

a = popt[0]
b = popt[1]
c = popt[2]
print a
print b
print c
