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
fileName = "nShellsEneOrbits";
fileNameSRG = "N26_srg";
fileNameFCIQ = "N26_fciqmc_data";
fileNameCCD = "nShellsEneAlpha0.3";

dataSRG = loadtxt(fileNameSRG)
dataFCIQ = loadtxt(fileNameFCIQ)
dataCCD = loadtxt(fileNameCCD)

orbitsSRG = 1/dataSRG[:,1]
orbitsFCIQ = 1/dataFCIQ[:,1]
orbitsCCD = 1/dataCCD[:,1]

eneSRG = dataSRG[:,3]
eneFCIQ = dataFCIQ[:,2]
eneCCD = dataCCD[:,3] 

dashes1 = [5,2,2,2];
dashes2 = [5,2,1,2,1,2];

spring_green = (0, 1.0, 0)

#     Plot
p1 = plot(orbitsSRG, eneSRG, color="red", linewidth=1.2, 
          linestyle=":", marker=".", label="SRG")
p1 = plot(orbitsFCIQ, eneFCIQ, color=spring_green, 
          linestyle="noline", marker=".", label="FCIQMC")
p3 = plot(orbitsCCD, eneCCD, color="blue", linewidth=1.2,
          linestyle="-.", marker=".", label="CCD")

plot([0, 0.021], [-0.217, -0.217], 'k', lw=1.2, linestyle="--",
     label="VMC, TDL")
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

xticks([0.0, 1./50, 1./98, 1./162, 1./394],
       [r'$\infty^{-1}$', r'$50^{-1}$', r'$98^{-1}$', r'$162^{-1}$', r'$394^{-1}$'])
#yticks([-0.217], [-0.217])

#     Limits on axes
#xlim(nShellsCCD.min(), nShellsCCD.max())
#xlim(0.0, orbitsSRG.max())
xlim(0.0, 0.021)
#ylim(ene_hf.min(), ene_hf.max())
#ylim(-0.23, -0.17)

#     Set axis labels
xlabel('$n_{\mathrm{basis}}^{-1}$', fontsize=16)
ylabel('Correlation energy per particle (Ry)', fontsize=16)
title("2D electron gas, $N = 26$, $r_{s} = 1.0$")

fig.tight_layout()
frame = axis.get_frame()
frame.set_linewidth(0.1)

show()

fig.savefig(fileName+'.pdf', format='pdf')
fig.savefig(fileName+'.eps', format='eps')
