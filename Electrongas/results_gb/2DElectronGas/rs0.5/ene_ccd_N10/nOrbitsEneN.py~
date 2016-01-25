from numpy import *
from matplotlib.pyplot import *
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#     Use Latex font as default
rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 14})
rc('legend',**{'fontsize':14})   

#     Initialize a figure
fig = figure(figsize=(6,4), dpi=60)

#     Read data file
fileName = "nOrbitsEneN10";
fileNameSRG10 = "N10_srg";
#fileNameSRG26 = "N26_srg";
#fileNameSRG42 = "N42_srg";
#fileNameSRG74 = "N74_srg";
#fileNameSRG98 = "N98_srg";
fileNameFCIQ = "N10_fciqmc";
fileNameCCD10 = "N10_ccd_alpha0.3"
# fileNameCCD74 = "../ene_ccd_N74/N74_ccd_alpha0.2";
# fileNameCCD114 = "../ene_ccd_N114/N114_ccd_alpha0.1";
# fileNameCCD138 = "../ene_ccd_N138/N138_ccd_alpha0.3";
# fileNameCCD178 = "../ene_ccd_N178/N178_ccd_alpha0.2";

dataSRG10 = loadtxt(fileNameSRG10)
#dataSRG74 = loadtxt(fileNameSRG74)
#dataSRG98 = loadtxt(fileNameSRG98)
dataFCIQ = loadtxt(fileNameFCIQ)
dataCCD10 = loadtxt(fileNameCCD10)
# dataCCD74 = loadtxt(fileNameCCD74)
# dataCCD114 = loadtxt(fileNameCCD114)
# dataCCD138 = loadtxt(fileNameCCD138)
# dataCCD178 = loadtxt(fileNameCCD178)

orbitsSRG10 = 1/dataSRG10[:,1]
#orbitsSRG74 = 1/dataSRG74[:,1]
#orbitsSRG98 = 1/dataSRG98[:,1]
orbitsFCIQ = 1/dataFCIQ[:,1]
orbitsCCD10 = 1/dataCCD10[:,1]
# orbitsCCD74 = 1/dataCCD74[:,1]
# orbitsCCD114 = 1/dataCCD114[:,1]
# orbitsCCD138 = 1/dataCCD138[:,1]
# orbitsCCD178 = 1/dataCCD178[:,1]

eneSRG10 = dataSRG10[:,3]
#eneSRG74 = dataSRG74[:,3]
#eneSRG98 = dataSRG98[:,3]
eneFCIQ = dataFCIQ[:,3]
eneCCD10 = dataCCD10[:,3]
# eneCCD74 = dataCCD74[:,3] 
# eneCCD114 = dataCCD114[:,3]
# eneCCD138 = dataCCD138[:,3]
# eneCCD178 = dataCCD178[:,3]

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
p1 = plot(orbitsSRG10, eneSRG10, color="red", linewidth=1.2, 
          linestyle="-", marker=".", label="SRG")
p4 = plot(orbitsCCD10, eneCCD10, color="blue", linewidth=1.2,
          linestyle="-", marker=".", label="CCD")
#p2 = plot(orbitsSRG98, eneSRG98, color="red", linewidth=1.2, 
#          linestyle="-", marker=".", label="SRG, $N = 98$")
p3 = plot(orbitsFCIQ, eneFCIQ, color=spring_green, 
          linestyle="noline", marker=".", label="FCIQMC")
# p5 = plot(orbitsCCD114, eneCCD114, color="blue", linewidth=1.2,
#           linestyle="-.", marker=".", label="CCD, $N = 114$")
# p6 = plot(orbitsCCD138, eneCCD138, color="blue", linewidth=1.2,
#           linestyle="--", marker=".", label="CCD, $N = 138$")
# p7 = plot(orbitsCCD178, eneCCD178, color="blue", linewidth=1.2,
#           linestyle="-", marker=".", label="CCD, $N = 178$")

#plot([0, 0.005], [-0.217, -0.217], 'k', lw=1.2, linestyle="--",
#     label="VMC, TDL")
# plot([0, N_ele.max()], [0.242, 0.242], 'k', lw=1.2, linestyle="-.")
# plot([0, N_ele.max()], [0.0945, 0.0945], 'k', lw=1.2, linestyle="-")


#     Legend
leg = legend(loc='upper left', labelspacing=0.1)
fr_leg = leg.get_frame()
fr_leg.set_lw(0.6)

#     Set minor and major ticks
#majorLocatorX = MultipleLocator(50)
#minorLocatorX = MultipleLocator(10)
majorLocatorY = MultipleLocator(0.02)
minorLocatorY = MultipleLocator(0.01)
axis = fig.gca()
#axis.xaxis.set_major_locator(majorLocatorX)
#axis.xaxis.set_minor_locator(minorLocatorX)
axis.yaxis.set_major_locator(majorLocatorY)
axis.yaxis.set_minor_locator(minorLocatorY)

xticks([0.0, 1./242, 1./98, 1./50, 1./26],
       [r'$\infty^{-1}$', r'$242^{-1}$', r'$98^{-1}$', 
        r'$50^{-1}$', r'$26^{-1}$'])
#yticks([-0.217], [-0.217])

#     Limits on axes
#xlim(nShellsCCD.min(), nShellsCCD.max())
#xlim(0.0, orbitsSRG.max())
#xlim(0.0, 0.0043)
#ylim(ene_hf.min(), ene_hf.max())
ylim(-0.21, -0.06)

#     Set axis labels
xlabel('Orbits$^{-1}$', fontsize=16)
ylabel('Correlation energy per particle (Ry)', fontsize=16)
title("2D electron gas, $r_{s} = 2.0$, $N = 10$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

#show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
